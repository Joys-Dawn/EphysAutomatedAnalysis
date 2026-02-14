"""
R bridge: call lmerTest + emmeans via subprocess for LMM inference.

All R communication flows through this module.  Each function writes a
DataFrame to a temp CSV, generates an R script, runs it via Rscript.exe,
and parses structured stdout lines back into Python dataclasses.
"""

import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

if sys.platform == 'win32':
    import winreg

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class OmnibusEffect:
    """One row from lmerTest::anova (Satterthwaite)."""
    effect: str
    f_value: float
    df1: float
    df2: float
    p_value: float


@dataclass
class PairwiseResult:
    """One row from emmeans::pairs."""
    contrast: str  # e.g. "A - B"
    level_a: str
    level_b: str
    estimate: float
    se: float
    df: float
    t_value: float
    p_value: float


@dataclass
class SimpleEffectResult:
    """One row from emmeans simple-effects analysis."""
    by_factor: str  # name of the conditioning factor (e.g. "factor1")
    by_level: str   # level of the conditioning factor
    contrast: str
    level_a: str
    level_b: str
    estimate: float
    se: float
    df: float
    t_value: float
    p_value: float


@dataclass
class EMMean:
    """One estimated marginal mean from emmeans."""
    factor: str   # "factor1" or "factor2"
    level: str
    emmean: float
    se: float
    df: float


@dataclass
class VarianceComponent:
    group: str
    variance: float
    std_dev: float


@dataclass
class LMMResult:
    """Full result of an R LMM fit."""
    omnibus: List[OmnibusEffect] = field(default_factory=list)
    pairwise: List[PairwiseResult] = field(default_factory=list)
    simple_effects: List[SimpleEffectResult] = field(default_factory=list)
    marginals: Dict[str, List[PairwiseResult]] = field(default_factory=dict)  # factor_name -> pairwise
    emmeans: List[EMMean] = field(default_factory=list)
    variance_components: List[VarianceComponent] = field(default_factory=list)
    singular: bool = False
    converged: bool = True


# ---------------------------------------------------------------------------
# R auto-detection (cached)
# ---------------------------------------------------------------------------

_RSCRIPT_PATH: Optional[str] = None
_R_CHECKED = False
_R_PACKAGES_OK: Optional[bool] = None


def find_r() -> Optional[str]:
    """Find Rscript on this system (Windows, macOS, Linux).  Cached after first call."""
    global _RSCRIPT_PATH, _R_CHECKED
    if _R_CHECKED:
        return _RSCRIPT_PATH
    _R_CHECKED = True

    # 1. R_HOME environment variable
    r_home = os.environ.get('R_HOME')
    if r_home:
        p = _rscript_from_home(r_home)
        if p:
            _RSCRIPT_PATH = p
            return _RSCRIPT_PATH

    # 2. Windows registry
    if sys.platform == 'win32':
        for hive in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
            for sub in (r'SOFTWARE\R-core\R64', r'SOFTWARE\R-core\R'):
                try:
                    key = winreg.OpenKey(hive, sub)
                    install_path, _ = winreg.QueryValueEx(key, 'InstallPath')
                    winreg.CloseKey(key)
                    p = _rscript_from_home(install_path)
                    if p:
                        _RSCRIPT_PATH = p
                        return _RSCRIPT_PATH
                except OSError:
                    continue

    # 3. PATH lookup (works on all platforms)
    p = shutil.which('Rscript')
    if p:
        _RSCRIPT_PATH = p
        return _RSCRIPT_PATH

    # 4. Common installation directories (platform-specific)
    if sys.platform == 'win32':
        search_patterns = (
            r'C:\R\R-*', r'C:\Program Files\R\R-*',
            os.path.expandvars(r'%LOCALAPPDATA%\Programs\R\R-*'),
        )
    elif sys.platform == 'darwin':
        search_patterns = (
            '/Library/Frameworks/R.framework/Versions/*',
            '/usr/local/Cellar/r/*/lib/R',
        )
    else:  # Linux
        search_patterns = (
            '/usr/lib/R',
            '/usr/local/lib/R',
            '/opt/R/*/lib/R',
        )

    for pattern in search_patterns:
        parent = Path(pattern).parent
        if parent.exists():
            dirs = sorted(parent.glob(Path(pattern).name), reverse=True)
            for d in dirs:
                p = _rscript_from_home(str(d))
                if p:
                    _RSCRIPT_PATH = p
                    return _RSCRIPT_PATH

    logger.warning("R installation not found. LMM will fall back to classical tests.")
    return None


def _rscript_from_home(r_home: str) -> Optional[str]:
    """Return path to Rscript inside an R home directory, or None."""
    for sub in ('bin/x64/Rscript.exe', 'bin/Rscript.exe',
                'bin/x64/Rscript', 'bin/Rscript'):
        p = os.path.join(r_home, sub)
        if os.path.isfile(p):
            return p
    return None


_REQUIRED_R_PACKAGES = ("lmerTest", "emmeans")


def _r_packages_present(rscript: str) -> bool:
    """Return True if all required R packages are importable."""
    check = " && ".join(
        f'requireNamespace("{pkg}", quietly=TRUE)' for pkg in _REQUIRED_R_PACKAGES
    )
    try:
        res = subprocess.run(
            [rscript, '-e', f'cat({check})'],
            capture_output=True, text=True, timeout=15,
        )
        return res.stdout.strip() == 'TRUE'
    except Exception:
        return False


def _install_r_packages(rscript: str) -> bool:
    """Install required R packages to user library. Returns True on success."""
    pkgs = ", ".join(f'"{pkg}"' for pkg in _REQUIRED_R_PACKAGES)
    # Install to user library so no admin/root is needed
    script = (
        'lib <- Sys.getenv("R_LIBS_USER"); '
        'if (lib == "") lib <- file.path(Sys.getenv("HOME"), "R", "library"); '
        'dir.create(lib, recursive=TRUE, showWarnings=FALSE); '
        f'install.packages(c({pkgs}), lib=lib, repos="https://cloud.r-project.org", quiet=TRUE); '
        f'cat(all(sapply(c({pkgs}), requireNamespace, quietly=TRUE)))'
    )
    try:
        logger.info("Installing required R packages (lmerTest, emmeans)...")
        res = subprocess.run(
            [rscript, '-e', script],
            capture_output=True, text=True, timeout=300,
        )
        if res.stdout.strip().endswith('TRUE'):
            logger.info("R packages installed successfully.")
            return True
        logger.warning(f"R package installation did not succeed. stderr: {res.stderr[:500]}")
        return False
    except subprocess.TimeoutExpired:
        logger.warning("R package installation timed out.")
        return False
    except Exception as e:
        logger.warning(f"R package installation failed: {e}")
        return False


def check_r_packages() -> bool:
    """Verify lmerTest and emmeans are installed, auto-installing if needed."""
    global _R_PACKAGES_OK
    if _R_PACKAGES_OK is not None:
        return _R_PACKAGES_OK

    rscript = find_r()
    if rscript is None:
        _R_PACKAGES_OK = False
        return False

    if _r_packages_present(rscript):
        _R_PACKAGES_OK = True
        return True

    # Attempt auto-install
    if _install_r_packages(rscript) and _r_packages_present(rscript):
        _R_PACKAGES_OK = True
        return True

    _R_PACKAGES_OK = False
    logger.warning("R packages lmerTest/emmeans not available. LMM disabled.")
    return False


def is_r_available() -> bool:
    """One-shot check: R installed AND required packages present."""
    return find_r() is not None and check_r_packages()


# ---------------------------------------------------------------------------
# Core subprocess runner
# ---------------------------------------------------------------------------

def _run_r_script(script: str, timeout: int = 60) -> Optional[str]:
    """Write *script* to a temp file, run via Rscript, return stdout or None."""
    rscript = find_r()
    if rscript is None:
        return None

    fd, path = tempfile.mkstemp(suffix='.R', prefix='lmm_')
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(script)
        result = subprocess.run(
            [rscript, '--vanilla', path],
            capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode != 0:
            # Filter out harmless loading messages
            stderr = '\n'.join(
                l for l in result.stderr.splitlines()
                if not any(x in l for x in ('Loading', 'Attaching', 'package', 'mask'))
            ).strip()
            if stderr:
                logger.debug(f"R stderr: {stderr[:300]}")
            return None
        return result.stdout
    except subprocess.TimeoutExpired:
        logger.warning("R script timed out")
        return None
    except Exception as e:
        logger.warning(f"R subprocess error: {e}")
        return None
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


def _write_temp_csv(data: pd.DataFrame) -> str:
    """Write DataFrame to temp CSV and return the path (forward-slashes for R)."""
    fd, path = tempfile.mkstemp(suffix='.csv', prefix='lmm_data_')
    os.close(fd)
    data.to_csv(path, index=False)
    return path.replace('\\', '/')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fit_lmm_anova(
    data: pd.DataFrame,
    formula: str,
    random_effects: str,
    ddf: str = "Satterthwaite",
) -> Optional[LMMResult]:
    """Fit an LMM and return Satterthwaite omnibus F-tests.

    Parameters
    ----------
    data : DataFrame with columns referenced in formula/random_effects.
    formula : R fixed-effects formula, e.g. ``"y ~ Group"``
    random_effects : R random-effects terms, e.g. ``"(1|Mouse_ID)"`` or
                     ``"(1|Mouse_ID) + (1|Cell_ID)"``
    ddf : degrees-of-freedom method (default ``"Satterthwaite"``).

    Returns LMMResult with omnibus populated, or None on failure/singular.
    """
    csv_path = _write_temp_csv(data)
    try:
        full_formula = f"{formula} + {random_effects}"
        script = _ANOVA_TEMPLATE.format(csv=csv_path, formula=full_formula, ddf=ddf)
        stdout = _run_r_script(script)
        if stdout is None:
            return None
        return _parse_anova_output(stdout)
    finally:
        _cleanup_csv(csv_path)


def fit_lmm_pairwise(
    data: pd.DataFrame,
    formula: str,
    random_effects: str,
    factor: str,
    ddf: str = "Satterthwaite",
) -> Optional[LMMResult]:
    """Fit LMM, run emmeans pairwise on *factor*, return omnibus + pairwise."""
    csv_path = _write_temp_csv(data)
    try:
        full_formula = f"{formula} + {random_effects}"
        script = _PAIRWISE_TEMPLATE.format(
            csv=csv_path, formula=full_formula, factor=factor, ddf=ddf,
        )
        stdout = _run_r_script(script)
        if stdout is None:
            return None
        return _parse_pairwise_output(stdout)
    finally:
        _cleanup_csv(csv_path)


def fit_lmm_factorial(
    data: pd.DataFrame,
    formula: str,
    random_effects: str,
    factor1: str,
    factor2: str,
    ddf: str = "Satterthwaite",
) -> Optional[LMMResult]:
    """Fit factorial LMM, return omnibus (3 effects) + simple effects."""
    csv_path = _write_temp_csv(data)
    try:
        full_formula = f"{formula} + {random_effects}"
        script = _FACTORIAL_TEMPLATE.format(
            csv=csv_path, formula=full_formula,
            factor1=factor1, factor2=factor2, ddf=ddf,
        )
        stdout = _run_r_script(script)
        if stdout is None:
            return None
        return _parse_factorial_output(stdout)
    finally:
        _cleanup_csv(csv_path)




# ---------------------------------------------------------------------------
# R script templates
# ---------------------------------------------------------------------------

_FIT_BLOCK = '''\
conv_warn <- NULL
fit <- tryCatch(
    withCallingHandlers(
        lmer({formula}, data=df, REML=TRUE),
        warning = function(w) {{
            if (grepl("converge|gradient|singular", conditionMessage(w), ignore.case=TRUE)) {{
                conv_warn <<- conditionMessage(w)
            }}
            invokeRestart("muffleWarning")
        }}
    ),
    error = function(e) NULL
)

if (is.null(fit)) {{
    cat("STATUS: FIT_FAILED\\n")
}} else {{
    is_sing <- isSingular(fit)
    conv_msgs <- fit@optinfo$conv$lme4$messages
    has_conv_issue <- !is.null(conv_warn) || !is.null(conv_msgs)
    if (is_sing) {{
        cat("STATUS: SINGULAR\\n")
    }} else if (has_conv_issue) {{
        cat("STATUS: CONVERGENCE_WARNING\\n")
    }} else {{
        cat("STATUS: OK\\n")
    }}
'''

_ANOVA_TEMPLATE = '''\
suppressPackageStartupMessages({{
    library(lmerTest)
}})

df <- read.csv("{csv}")
for (col in names(df)) {{
    if (is.character(df[[col]])) df[[col]] <- factor(df[[col]])
}}

''' + _FIT_BLOCK + '''
    # Variance components
    vc <- as.data.frame(VarCorr(fit))
    for (i in 1:nrow(vc)) {{
        cat(sprintf("VC: group=%s variance=%.10e stddev=%.10e\\n",
            vc$grp[i], vc$vcov[i], vc$sdcor[i]))
    }}

    # Satterthwaite ANOVA
    aov <- anova(fit, ddf="{ddf}")
    for (i in 1:nrow(aov)) {{
        cat(sprintf("OMNIBUS: effect=%s F=%.10e df1=%.4f df2=%.4f p=%.10e\\n",
            rownames(aov)[i], aov[["F value"]][i],
            aov[["NumDF"]][i], aov[["DenDF"]][i], aov[["Pr(>F)"]][i]))
    }}
}}
'''

_PAIRWISE_TEMPLATE = '''\
suppressPackageStartupMessages({{
    library(lmerTest)
    library(emmeans)
}})

df <- read.csv("{csv}")
for (col in names(df)) {{
    if (is.character(df[[col]])) df[[col]] <- factor(df[[col]])
}}

''' + _FIT_BLOCK + '''
    # Variance components
    vc <- as.data.frame(VarCorr(fit))
    for (i in 1:nrow(vc)) {{
        cat(sprintf("VC: group=%s variance=%.10e stddev=%.10e\\n",
            vc$grp[i], vc$vcov[i], vc$sdcor[i]))
    }}

    # Satterthwaite ANOVA
    aov <- anova(fit, ddf="{ddf}")
    for (i in 1:nrow(aov)) {{
        cat(sprintf("OMNIBUS: effect=%s F=%.10e df1=%.4f df2=%.4f p=%.10e\\n",
            rownames(aov)[i], aov[["F value"]][i],
            aov[["NumDF"]][i], aov[["DenDF"]][i], aov[["Pr(>F)"]][i]))
    }}

    # emmeans pairwise
    emm <- emmeans(fit, as.formula(paste("~", "{factor}")), lmer.df="{ddf}")
    pw <- as.data.frame(pairs(emm, adjust="none"))
    for (i in 1:nrow(pw)) {{
        cat(sprintf("PAIR: contrast=%s estimate=%.10e SE=%.10e df=%.4f t=%.10e p=%.10e\\n",
            pw$contrast[i], pw$estimate[i], pw$SE[i],
            pw$df[i], pw$t.ratio[i], pw$p.value[i]))
    }}
}}
'''

_FACTORIAL_TEMPLATE = '''\
suppressPackageStartupMessages({{
    library(lmerTest)
    library(emmeans)
}})

df <- read.csv("{csv}")
for (col in names(df)) {{
    if (is.character(df[[col]])) df[[col]] <- factor(df[[col]])
}}

''' + _FIT_BLOCK + '''
    # Variance components
    vc <- as.data.frame(VarCorr(fit))
    for (i in 1:nrow(vc)) {{
        cat(sprintf("VC: group=%s variance=%.10e stddev=%.10e\\n",
            vc$grp[i], vc$vcov[i], vc$sdcor[i]))
    }}

    # Satterthwaite ANOVA (all effects)
    aov <- anova(fit, ddf="{ddf}")
    for (i in 1:nrow(aov)) {{
        cat(sprintf("OMNIBUS: effect=%s F=%.10e df1=%.4f df2=%.4f p=%.10e\\n",
            rownames(aov)[i], aov[["F value"]][i],
            aov[["NumDF"]][i], aov[["DenDF"]][i], aov[["Pr(>F)"]][i]))
    }}

    # Simple effects: compare factor1 at each level of factor2
    emm <- emmeans(fit, as.formula(paste("~", "{factor1}", "|", "{factor2}")),
                   lmer.df="{ddf}")
    se_df <- as.data.frame(pairs(emm, adjust="none"))
    for (i in 1:nrow(se_df)) {{
        cat(sprintf("SIMPLE: by_factor=factor2 by=%s contrast=%s estimate=%.10e SE=%.10e df=%.4f t=%.10e p=%.10e\\n",
            se_df${factor2}[i], se_df$contrast[i], se_df$estimate[i], se_df$SE[i],
            se_df$df[i], se_df$t.ratio[i], se_df$p.value[i]))
    }}

    # Simple effects: compare factor2 at each level of factor1
    emm2 <- emmeans(fit, as.formula(paste("~", "{factor2}", "|", "{factor1}")),
                    lmer.df="{ddf}")
    se_df2 <- as.data.frame(pairs(emm2, adjust="none"))
    for (i in 1:nrow(se_df2)) {{
        cat(sprintf("SIMPLE: by_factor=factor1 by=%s contrast=%s estimate=%.10e SE=%.10e df=%.4f t=%.10e p=%.10e\\n",
            se_df2${factor1}[i], se_df2$contrast[i], se_df2$estimate[i], se_df2$SE[i],
            se_df2$df[i], se_df2$t.ratio[i], se_df2$p.value[i]))
    }}

    # Marginal comparisons: factor1 averaged over factor2
    emm_m1 <- emmeans(fit, as.formula(paste("~", "{factor1}")), lmer.df="{ddf}")
    emm_m1_df <- as.data.frame(emm_m1)
    for (i in 1:nrow(emm_m1_df)) {{
        cat(sprintf("EMM: factor=factor1 level=%s emmean=%.10e SE=%.10e df=%.4f\\n",
            emm_m1_df[i, 1], emm_m1_df$emmean[i], emm_m1_df$SE[i], emm_m1_df$df[i]))
    }}
    pw_m1 <- as.data.frame(pairs(emm_m1, adjust="none"))
    for (i in 1:nrow(pw_m1)) {{
        cat(sprintf("MARGINAL: factor=factor1 contrast=%s estimate=%.10e SE=%.10e df=%.4f t=%.10e p=%.10e\\n",
            pw_m1$contrast[i], pw_m1$estimate[i], pw_m1$SE[i],
            pw_m1$df[i], pw_m1$t.ratio[i], pw_m1$p.value[i]))
    }}

    # Marginal comparisons: factor2 averaged over factor1
    emm_m2 <- emmeans(fit, as.formula(paste("~", "{factor2}")), lmer.df="{ddf}")
    emm_m2_df <- as.data.frame(emm_m2)
    for (i in 1:nrow(emm_m2_df)) {{
        cat(sprintf("EMM: factor=factor2 level=%s emmean=%.10e SE=%.10e df=%.4f\\n",
            emm_m2_df[i, 1], emm_m2_df$emmean[i], emm_m2_df$SE[i], emm_m2_df$df[i]))
    }}
    pw_m2 <- as.data.frame(pairs(emm_m2, adjust="none"))
    for (i in 1:nrow(pw_m2)) {{
        cat(sprintf("MARGINAL: factor=factor2 contrast=%s estimate=%.10e SE=%.10e df=%.4f t=%.10e p=%.10e\\n",
            pw_m2$contrast[i], pw_m2$estimate[i], pw_m2$SE[i],
            pw_m2$df[i], pw_m2$t.ratio[i], pw_m2$p.value[i]))
    }}
}}
'''



# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------

def _cleanup_csv(path: str):
    """Remove a temp CSV file."""
    try:
        os.unlink(path)
    except OSError:
        pass


def _parse_status(stdout: str) -> tuple:
    """Return (ok: bool, singular: bool) from STATUS line."""
    for line in stdout.splitlines():
        if line.startswith('STATUS:'):
            status = line.split(':', 1)[1].strip()
            if status == 'FIT_FAILED':
                return False, False
            if status == 'SINGULAR':
                return False, True  # singular -> fallback
            if status == 'CONVERGENCE_WARNING':
                logger.warning("LMM convergence warning -- falling back to classical")
                return False, False  # convergence issue -> fallback
            return True, False
    return False, False


def _parse_vc(stdout: str) -> List[VarianceComponent]:
    """Parse VC: lines into VarianceComponent list."""
    result = []
    for line in stdout.splitlines():
        if line.startswith('VC:'):
            parts = _parse_kv(line[3:].strip())
            result.append(VarianceComponent(
                group=parts.get('group', ''),
                variance=float(parts.get('variance', 0)),
                std_dev=float(parts.get('stddev', 0)),
            ))
    return result


def _parse_omnibus(stdout: str) -> List[OmnibusEffect]:
    """Parse OMNIBUS: lines."""
    result = []
    for line in stdout.splitlines():
        if line.startswith('OMNIBUS:'):
            parts = _parse_kv(line[8:].strip())
            result.append(OmnibusEffect(
                effect=parts.get('effect', ''),
                f_value=float(parts.get('F', 0)),
                df1=float(parts.get('df1', 0)),
                df2=float(parts.get('df2', 0)),
                p_value=float(parts.get('p', 1)),
            ))
    return result


def _parse_pairs(stdout: str) -> List[PairwiseResult]:
    """Parse PAIR: lines into PairwiseResult list."""
    result = []
    for line in stdout.splitlines():
        if line.startswith('PAIR:'):
            parts = _parse_kv(line[5:].strip())
            contrast = parts.get('contrast', '')
            la, lb = _split_contrast(contrast)
            result.append(PairwiseResult(
                contrast=contrast,
                level_a=la,
                level_b=lb,
                estimate=float(parts.get('estimate', 0)),
                se=float(parts.get('SE', 0)),
                df=float(parts.get('df', 0)),
                t_value=float(parts.get('t', 0)),
                p_value=float(parts.get('p', 1)),
            ))
    return result


def _parse_simple_effects(stdout: str) -> List[SimpleEffectResult]:
    """Parse SIMPLE: lines."""
    result = []
    for line in stdout.splitlines():
        if line.startswith('SIMPLE:'):
            parts = _parse_kv(line[7:].strip())
            contrast = parts.get('contrast', '')
            la, lb = _split_contrast(contrast)
            result.append(SimpleEffectResult(
                by_factor=parts.get('by_factor', ''),
                by_level=parts.get('by', ''),
                contrast=contrast,
                level_a=la,
                level_b=lb,
                estimate=float(parts.get('estimate', 0)),
                se=float(parts.get('SE', 0)),
                df=float(parts.get('df', 0)),
                t_value=float(parts.get('t', 0)),
                p_value=float(parts.get('p', 1)),
            ))
    return result


def _parse_marginals(stdout: str) -> Dict[str, List[PairwiseResult]]:
    """Parse MARGINAL: lines into {factor_name: [PairwiseResult, ...]}."""
    result: Dict[str, List[PairwiseResult]] = {}
    for line in stdout.splitlines():
        if line.startswith('MARGINAL:'):
            parts = _parse_kv(line[9:].strip())
            factor = parts.get('factor', '')
            contrast = parts.get('contrast', '')
            la, lb = _split_contrast(contrast)
            pw = PairwiseResult(
                contrast=contrast,
                level_a=la,
                level_b=lb,
                estimate=float(parts.get('estimate', 0)),
                se=float(parts.get('SE', 0)),
                df=float(parts.get('df', 0)),
                t_value=float(parts.get('t', 0)),
                p_value=float(parts.get('p', 1)),
            )
            result.setdefault(factor, []).append(pw)
    return result


def _parse_emmeans(stdout: str) -> List[EMMean]:
    """Parse EMM: lines into EMMean list."""
    result = []
    for line in stdout.splitlines():
        if line.startswith('EMM:'):
            parts = _parse_kv(line[4:].strip())
            result.append(EMMean(
                factor=parts.get('factor', ''),
                level=parts.get('level', ''),
                emmean=float(parts.get('emmean', 0)),
                se=float(parts.get('SE', 0)),
                df=float(parts.get('df', 0)),
            ))
    return result


def _parse_anova_output(stdout: str) -> Optional[LMMResult]:
    """Parse output from the anova-only template."""
    ok, singular = _parse_status(stdout)
    if not ok:
        return None
    return LMMResult(
        omnibus=_parse_omnibus(stdout),
        variance_components=_parse_vc(stdout),
        singular=singular,
    )


def _parse_pairwise_output(stdout: str) -> Optional[LMMResult]:
    """Parse output from the pairwise template."""
    ok, singular = _parse_status(stdout)
    if not ok:
        return None
    return LMMResult(
        omnibus=_parse_omnibus(stdout),
        pairwise=_parse_pairs(stdout),
        variance_components=_parse_vc(stdout),
        singular=singular,
    )


def _parse_factorial_output(stdout: str) -> Optional[LMMResult]:
    """Parse output from the factorial template (omnibus + simple + marginals + emmeans)."""
    ok, singular = _parse_status(stdout)
    if not ok:
        return None
    return LMMResult(
        omnibus=_parse_omnibus(stdout),
        simple_effects=_parse_simple_effects(stdout),
        marginals=_parse_marginals(stdout),
        emmeans=_parse_emmeans(stdout),
        variance_components=_parse_vc(stdout),
        singular=singular,
    )


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _parse_kv(text: str) -> Dict[str, str]:
    """Parse ``key=value key2=value2 ...`` into a dict.

    Handles values containing spaces (e.g. contrast names like
    ``"32 Scn1a - 37 Scn1a"``) by scanning for ``word=`` boundaries.
    """
    result = {}
    # Find all key=value pairs.  Keys are single words; values run until
    # the next ``word=`` boundary or end of string.
    pattern = re.compile(r'(\w+)=(.*?)(?=\s+\w+=|$)')
    for m in pattern.finditer(text):
        result[m.group(1)] = m.group(2).strip()
    return result


def _split_contrast(contrast: str) -> tuple:
    """Split ``"A - B"`` into ``("A", "B")``.  Handles multi-word level names."""
    parts = contrast.split(' - ', 1)
    if len(parts) == 2:
        return parts[0].strip(), parts[1].strip()
    return contrast.strip(), ''
