"""
LMM-enhanced mixed ANOVA for between x within factorial designs.

Subclasses MixedANOVA: fits an LMM with Mouse_ID and Cell_ID as nested
random intercepts via R's lmerTest (Satterthwaite df) for omnibus effects,
and emmeans for simple effects and marginal comparisons.  Falls back to
super() (classical pg.mixed_anova) when R is unavailable or the model
fails/is singular.
"""

import os
import math
import logging
import pandas as pd
from typing import Optional, List, Dict

from ..statistical_analysis.tests.mixed_anova import MixedANOVA
from ..shared.data_models import ExperimentalDesign, StatisticalResult
from ..shared.utils import categorize_measurement
from .mouse_log import load_mouse_log
from . import r_bridge

logger = logging.getLogger(__name__)


class LMMMixedANOVA(MixedANOVA):
    """Mixed ANOVA with LMM enhancement for mouse clustering."""

    def __init__(self, mouse_log_path: str):
        super().__init__()
        self.name = "Linear Mixed Model (mixed)"
        self.mouse_log = load_mouse_log(mouse_log_path)
        # measurement -> pre-computed LMM simple effects
        self._lmm_simple_effects: Dict[str, List[StatisticalResult]] = {}
        # measurement -> full LMMResult from factorial fit
        self._lmm_results: Dict[str, 'r_bridge.LMMResult'] = {}
        # measurement -> sub DataFrame used for LMM fitting
        self._lmm_sub_data: Dict[str, pd.DataFrame] = {}

    # ------------------------------------------------------------------
    # Data enrichment: add Mouse_ID after parent builds unified_df
    # ------------------------------------------------------------------

    def _create_unified_dataframe(self, group_data, design, manifest):
        unified = super()._create_unified_dataframe(group_data, design, manifest)
        if 'filename' in unified.columns:
            unified['Mouse_ID'] = unified['filename'].apply(
                lambda f: self.mouse_log.get(os.path.basename(str(f).strip()))
            )
            missing = unified['Mouse_ID'].isna().sum()
            if missing > 0:
                logger.warning(f"Mixed data: {missing}/{len(unified)} files not in mouse log, dropping")
                unified = unified.dropna(subset=['Mouse_ID'])
        return unified

    # ------------------------------------------------------------------
    # Omnibus: try LMM via R for each effect, fall back to super()
    # ------------------------------------------------------------------

    def _run_single_mixed_model(self, column: str, unified_df: pd.DataFrame,
                                design: ExperimentalDesign) -> Optional[List[StatisticalResult]]:
        """Try nested LMM + Satterthwaite for 3 effects; fall back to classical."""

        if not r_bridge.is_r_available():
            return super()._run_single_mixed_model(column, unified_df, design)

        if 'Mouse_ID' not in unified_df.columns:
            return super()._run_single_mixed_model(column, unified_df, design)

        if column not in unified_df.columns:
            return None

        # Build working sub-DataFrame with LMM-friendly column names
        cols = ['Between_Factor', 'Within_Factor', 'Subject_ID', 'Mouse_ID', column]
        df_clean = unified_df[cols].copy().dropna(subset=[column])

        # Collapse to one observation per subject/condition
        if not df_clean.empty:
            df_clean = (
                df_clean.groupby(['Subject_ID', 'Between_Factor', 'Within_Factor', 'Mouse_ID'],
                                 as_index=False)[column].mean()
            )

        # Require at least 2 observations per cell and 2 mice
        cell_counts = df_clean.groupby(['Between_Factor', 'Within_Factor']).size()
        if cell_counts.min() < 2 or df_clean['Mouse_ID'].nunique() < 2:
            return super()._run_single_mixed_model(column, unified_df, design)

        # Rename for R
        sub = df_clean.rename(columns={
            'Between_Factor': 'Between',
            'Within_Factor': 'Within',
            'Subject_ID': 'Cell_ID',
            column: 'y',
        })
        sub['y'] = pd.to_numeric(sub['y'], errors='coerce')
        sub = sub.dropna(subset=['y'])

        # Ensure Cell_ID is unique across mice to avoid crossed RE misspecification
        sub['Cell_ID'] = sub['Mouse_ID'].astype(str) + '_' + sub['Cell_ID'].astype(str)

        between_name = design.between_factor_name or "Between"
        within_name = design.within_factor_name or "Within"

        # Fit factorial LMM via R with nested RE
        result = r_bridge.fit_lmm_factorial(
            sub, 'y ~ Between * Within', '(1|Mouse_ID) + (1|Cell_ID)',
            'Between', 'Within',
        )
        if result is None:
            logger.info(f"LMM failed for '{column}', falling back to classical mixed ANOVA")
            return super()._run_single_mixed_model(column, unified_df, design)

        if not result.omnibus:
            logger.warning(f"LMM returned no omnibus for '{column}', falling back to classical")
            return super()._run_single_mixed_model(column, unified_df, design)

        measurement_type = categorize_measurement(column)

        # Store data and result for post-hocs
        self._lmm_results[column] = result
        self._lmm_sub_data[column] = sub

        # Map R omnibus effect names to design labels
        effect_map = {}
        for eff in result.omnibus:
            if eff.effect == 'Between':
                effect_map['between'] = eff.p_value
            elif eff.effect == 'Within':
                effect_map['within'] = eff.p_value
            elif 'Between' in eff.effect and 'Within' in eff.effect:
                effect_map['interaction'] = eff.p_value

        # Build omnibus results matching classical format
        effects = []
        for p_key, test_name, g1_name in [
            ('between', f"{self.name} - Between Effect", between_name),
            ('within', f"{self.name} - Within Effect", within_name),
            ('interaction', f"{self.name} - Interaction", f"{between_name} × {within_name}"),
        ]:
            effects.append(StatisticalResult(
                test_name=test_name,
                measurement=column,
                group1_name=g1_name,
                group1_mean=0, group1_stderr=0, group1_n=0,
                group2_name="",
                group2_mean=0, group2_stderr=0, group2_n=0,
                p_value=effect_map.get(p_key, 1.0),
                measurement_type=measurement_type,
            ))

        # Pre-compute simple effects from R emmeans
        simple_effects = []
        for se in result.simple_effects:
            ctx = se.by_level
            la, lb = se.level_a, se.level_b

            # by_factor tells us which factor is held constant
            if se.by_factor == 'factor2':
                # factor2 = Within held constant, comparing Between levels
                d1 = sub.loc[(sub['Between'] == la) & (sub['Within'] == ctx), 'y']
                d2 = sub.loc[(sub['Between'] == lb) & (sub['Within'] == ctx), 'y']
                g1n = f"{ctx}: {la}"
                g2n = f"{ctx}: {lb}"
                test_type = "LMM contrast (between-subject simple effect)"
            else:
                # factor1 = Between held constant, comparing Within levels
                d1 = sub.loc[(sub['Between'] == ctx) & (sub['Within'] == la), 'y']
                d2 = sub.loc[(sub['Between'] == ctx) & (sub['Within'] == lb), 'y']
                g1n = f"{ctx}: {la}"
                g2n = f"{ctx}: {lb}"
                test_type = "LMM contrast (within-subject simple effect)"

            simple_effects.append(StatisticalResult(
                test_name=test_type,
                measurement=column,
                group1_name=g1n,
                group1_mean=d1.mean() if len(d1) > 0 else 0.0,
                group1_stderr=d1.std(ddof=1) / math.sqrt(len(d1)) if len(d1) > 1 else 0.0,
                group1_n=len(d1),
                group2_name=g2n,
                group2_mean=d2.mean() if len(d2) > 0 else 0.0,
                group2_stderr=d2.std(ddof=1) / math.sqrt(len(d2)) if len(d2) > 1 else 0.0,
                group2_n=len(d2),
                p_value=se.p_value,
                measurement_type=measurement_type,
            ))
        self._lmm_simple_effects[column] = simple_effects

        return effects

    # ------------------------------------------------------------------
    # Simple effects: return LMM emmeans results or delegate to super()
    # ------------------------------------------------------------------

    def _run_posthoc_for_measurement(self, measurement: str, unified_df: pd.DataFrame,
                                     design: ExperimentalDesign,
                                     sig_dict=None) -> List[StatisticalResult]:
        if measurement in self._lmm_simple_effects:
            return self._lmm_simple_effects[measurement]
        return super()._run_posthoc_for_measurement(measurement, unified_df, design, sig_dict)

    # ------------------------------------------------------------------
    # Marginal comparisons: LMM emmeans marginal or delegate to super()
    # ------------------------------------------------------------------

    def _run_marginal_posthoc(self, measurement: str, unified_df: pd.DataFrame,
                              design: ExperimentalDesign, factor_key: str) -> List[StatisticalResult]:

        lmm_result = self._lmm_results.get(measurement)
        sub = self._lmm_sub_data.get(measurement)
        if lmm_result is None or sub is None:
            return super()._run_marginal_posthoc(measurement, unified_df, design, factor_key)

        between_name = design.between_factor_name or "Between"
        within_name = design.within_factor_name or "Within"

        # R call used factor1='Between', factor2='Within'
        if factor_key == 'between':
            target_col = 'Between'
            factor_label = between_name
            marginal_key = 'factor1'
        else:
            target_col = 'Within'
            factor_label = within_name
            marginal_key = 'factor2'

        marginal_pairwise = lmm_result.marginals.get(marginal_key, [])
        if not marginal_pairwise:
            return super()._run_marginal_posthoc(measurement, unified_df, design, factor_key)

        measurement_type = categorize_measurement(measurement)

        # Build lookup for model-based emmeans (mean + SE) per level
        emm_lookup = {}
        for emm in lmm_result.emmeans:
            if emm.factor == marginal_key:
                emm_lookup[emm.level] = emm

        # Total N per level from raw data
        n_per_level = sub.groupby(target_col)['y'].count().to_dict()

        results = []
        for pw in marginal_pairwise:
            test_name = ("LMM marginal contrast (within-subject)"
                         if factor_key == 'within'
                         else "LMM marginal contrast (between-subject)")

            emm_a = emm_lookup.get(pw.level_a)
            emm_b = emm_lookup.get(pw.level_b)

            results.append(StatisticalResult(
                test_name=test_name,
                measurement=measurement,
                group1_name=f"{factor_label}={pw.level_a}",
                group1_mean=emm_a.emmean if emm_a else 0.0,
                group1_stderr=emm_a.se if emm_a else 0.0,
                group1_n=n_per_level.get(pw.level_a, 0),
                group2_name=f"{factor_label}={pw.level_b}",
                group2_mean=emm_b.emmean if emm_b else 0.0,
                group2_stderr=emm_b.se if emm_b else 0.0,
                group2_n=n_per_level.get(pw.level_b, 0),
                p_value=pw.p_value,
                measurement_type=measurement_type,
            ))

        return results
