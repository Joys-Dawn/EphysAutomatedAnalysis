"""
LMM-enhanced one-way ANOVA for 3+ independent groups.

Subclasses OneWayANOVA: fits an LMM with Mouse_ID as a random intercept
via R's lmerTest (Satterthwaite df) for omnibus, and emmeans for pairwise
post-hocs.  Falls back to super() (classical Welch ANOVA / Kruskal-Wallis)
when R is unavailable or the model fails/is singular.
"""

import os
import math
import logging
import pandas as pd
from typing import Optional, List, Dict

from ..statistical_analysis.tests.oneway_anova import OneWayANOVA
from ..shared.data_models import StatisticalResult
from ..shared.utils import categorize_measurement
from .mouse_log import load_mouse_log
from . import r_bridge

logger = logging.getLogger(__name__)


class LMMOneWayANOVA(OneWayANOVA):
    """One-way ANOVA with LMM enhancement for mouse clustering."""

    def __init__(self, mouse_log_path: str):
        super().__init__()
        self.name = "Linear Mixed Model (one-way)"
        self.mouse_log = load_mouse_log(mouse_log_path)
        # measurement -> list of StatisticalResult (pre-computed pairwise from LMM)
        self._lmm_pairwise: Dict[str, List[StatisticalResult]] = {}

    # ------------------------------------------------------------------
    # Data enrichment
    # ------------------------------------------------------------------

    def _load_combined_data(self, group_name: str, base_path: str) -> pd.DataFrame:
        """Load CSV data then enrich with Mouse_ID from the mouse log."""
        df = super()._load_combined_data(group_name, base_path)
        if df.empty or 'filename' not in df.columns:
            return df
        df['Mouse_ID'] = df['filename'].apply(
            lambda f: self.mouse_log.get(os.path.basename(str(f).strip()))
        )
        missing = df['Mouse_ID'].isna().sum()
        if missing > 0:
            logger.warning(
                f"Group '{group_name}': {missing}/{len(df)} files not in mouse log, dropping"
            )
            df = df.dropna(subset=['Mouse_ID'])
        return df

    # ------------------------------------------------------------------
    # Per-measurement omnibus: try LMM, fall back to super()
    # ------------------------------------------------------------------

    def _run_single_anova(self, column: str, group_data: Dict[str, pd.DataFrame],
                          groups) -> Optional[StatisticalResult]:
        """Try LMM for omnibus; store pairwise if successful; else classical."""

        if not r_bridge.is_r_available():
            return super()._run_single_anova(column, group_data, groups)

        if not all('Mouse_ID' in gdf.columns for gdf in group_data.values()):
            return super()._run_single_anova(column, group_data, groups)

        # Merge into long-format DataFrame
        frames = []
        for gname, gdf in group_data.items():
            sub = gdf[['Mouse_ID', column]].dropna(subset=[column]).copy()
            sub['Group'] = gname
            frames.append(sub)
        merged = pd.concat(frames, ignore_index=True)
        merged.rename(columns={column: 'y'}, inplace=True)

        group_counts = merged.groupby('Group')['y'].count()
        valid_groups = group_counts[group_counts > 1].index.tolist()
        if len(valid_groups) < 2 or merged['Mouse_ID'].nunique() < 2:
            return super()._run_single_anova(column, group_data, groups)
        merged = merged[merged['Group'].isin(valid_groups)]

        # Fit LMM via R with emmeans pairwise
        result = r_bridge.fit_lmm_pairwise(
            merged, 'y ~ Group', '(1|Mouse_ID)', 'Group',
        )
        if result is None:
            logger.info(f"LMM failed for '{column}', falling back to classical ANOVA")
            return super()._run_single_anova(column, group_data, groups)

        # Omnibus p-value
        if not result.omnibus:
            logger.warning(f"LMM returned no omnibus for '{column}', falling back to classical")
            return super()._run_single_anova(column, group_data, groups)
        omnibus_p = result.omnibus[0].p_value

        # Store pairwise results for post-hoc phase
        levels = sorted(merged['Group'].unique().tolist())
        gdata = {g: merged.loc[merged['Group'] == g, 'y'] for g in levels}
        measurement_type = categorize_measurement(column)

        pw_results = []
        for pw in result.pairwise:
            ga = gdata.get(pw.level_a, pd.Series(dtype=float))
            gb = gdata.get(pw.level_b, pd.Series(dtype=float))
            pw_results.append(StatisticalResult(
                test_name=f"Pairwise LMM ({pw.level_a} vs {pw.level_b})",
                measurement=column,
                group1_name=pw.level_a,
                group1_mean=ga.mean() if len(ga) > 0 else 0.0,
                group1_stderr=ga.std(ddof=1) / math.sqrt(len(ga)) if len(ga) > 1 else 0.0,
                group1_n=len(ga),
                group2_name=pw.level_b,
                group2_mean=gb.mean() if len(gb) > 0 else 0.0,
                group2_stderr=gb.std(ddof=1) / math.sqrt(len(gb)) if len(gb) > 1 else 0.0,
                group2_n=len(gb),
                p_value=pw.p_value,
                measurement_type=measurement_type,
            ))
        self._lmm_pairwise[column] = pw_results

        # Build omnibus result
        first, second = levels[0], levels[1]
        return StatisticalResult(
            test_name=f"{self.name} (LMM)",
            measurement=column,
            group1_name=first,
            group1_mean=gdata[first].mean() if first in gdata and len(gdata[first]) > 0 else 0.0,
            group1_stderr=(gdata[first].std(ddof=1) / math.sqrt(len(gdata[first]))
                           if first in gdata and len(gdata[first]) > 1 else 0.0),
            group1_n=len(gdata[first]) if first in gdata else 0,
            group2_name=second,
            group2_mean=gdata[second].mean() if second in gdata and len(gdata[second]) > 0 else 0.0,
            group2_stderr=(gdata[second].std(ddof=1) / math.sqrt(len(gdata[second]))
                           if second in gdata and len(gdata[second]) > 1 else 0.0),
            group2_n=len(gdata[second]) if second in gdata else 0,
            p_value=omnibus_p,
            measurement_type=measurement_type,
        )

    # ------------------------------------------------------------------
    # Pairwise: return stored LMM pairwise or delegate to super()
    # ------------------------------------------------------------------

    def _run_pairwise_comparisons_if_significant(
        self, anova_results: List[StatisticalResult],
        group_data: Dict[str, pd.DataFrame], groups
    ) -> List[StatisticalResult]:
        """Return LMM emmeans pairwise for LMM measurements, classical for the rest."""

        significant_measurements = {
            r.measurement for r in anova_results
            if (r.corrected_p if r.corrected_p is not None else r.p_value) < 0.05
        }

        lmm_pairwise = []
        classical_anova_results = []

        for r in anova_results:
            if r.measurement not in significant_measurements:
                continue
            if r.measurement in self._lmm_pairwise:
                lmm_pairwise.extend(self._lmm_pairwise[r.measurement])
            else:
                classical_anova_results.append(r)

        classical_pairwise = []
        if classical_anova_results:
            classical_pairwise = super()._run_pairwise_comparisons_if_significant(
                classical_anova_results, group_data, groups
            )

        return lmm_pairwise + classical_pairwise
