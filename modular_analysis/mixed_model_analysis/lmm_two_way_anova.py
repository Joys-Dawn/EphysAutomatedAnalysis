"""
LMM-enhanced two-way ANOVA for N x M factorial designs.

Subclasses TwoWayANOVA: fits an LMM with Mouse_ID as a random intercept
via R's lmerTest (Satterthwaite df) for omnibus (3 effects), and emmeans
for simple effects and marginal comparisons.  Falls back to super()
(classical OLS ANOVA) when R is unavailable or the model fails/is singular.
"""

import os
import math
import logging
import pandas as pd
from typing import Optional, List, Dict

from ..statistical_analysis.tests.two_way_anova import TwoWayANOVA
from ..shared.data_models import ExperimentalDesign, StatisticalResult
from ..shared.utils import categorize_measurement
from .mouse_log import load_mouse_log
from . import r_bridge

logger = logging.getLogger(__name__)


class LMMTwoWayANOVA(TwoWayANOVA):
    """Two-way ANOVA with LMM enhancement for mouse clustering."""

    def __init__(self, mouse_log_path: str):
        super().__init__()
        self.name = "Linear Mixed Model (two-way)"
        self.mouse_log = load_mouse_log(mouse_log_path)
        # measurement -> pre-computed simple effects from LMM
        self._lmm_simple_effects: Dict[str, List[StatisticalResult]] = {}
        # measurement -> r_bridge LMMResult (for marginal post-hocs)
        self._lmm_results: Dict[str, r_bridge.LMMResult] = {}
        # measurement -> sub DataFrame used for LMM fitting
        self._lmm_sub_data: Dict[str, pd.DataFrame] = {}

    # ------------------------------------------------------------------
    # Data enrichment
    # ------------------------------------------------------------------

    def _load_combined_data(self, group_name: str, base_path: str) -> pd.DataFrame:
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
    # Omnibus: try LMM for the three effects, fall back to super()
    # ------------------------------------------------------------------

    def _run_single_two_way_anova(self, column: str, group_data: Dict[str, pd.DataFrame],
                                  design: ExperimentalDesign) -> Optional[List[StatisticalResult]]:
        """Try LMM for 3 factorial effects; fall back to classical OLS."""

        if not r_bridge.is_r_available():
            return super()._run_single_two_way_anova(column, group_data, design)

        if not all('Mouse_ID' in gdf.columns for gdf in group_data.values()):
            return super()._run_single_two_way_anova(column, group_data, design)

        # Build merged frame with factor1, factor2, Mouse_ID
        frames = []
        for group in design.groups:
            if group.name not in group_data:
                continue
            gdf = group_data[group.name]
            if column not in gdf.columns:
                continue
            chunk = gdf[['Mouse_ID', column]].dropna(subset=[column]).copy()
            chunk.rename(columns={column: 'y'}, inplace=True)
            chunk['factor1'] = design.factor_mapping[group.name]['factor1']
            chunk['factor2'] = design.factor_mapping[group.name]['factor2']
            chunk['group'] = group.name
            frames.append(chunk)

        if not frames:
            return super()._run_single_two_way_anova(column, group_data, design)

        sub = pd.concat(frames, ignore_index=True)
        sub['y'] = pd.to_numeric(sub['y'], errors='coerce')
        sub = sub.dropna(subset=['y'])

        # Require at least 2 observations per cell and 2 mice
        cell_counts = sub.groupby(['factor1', 'factor2']).size()
        if cell_counts.min() < 2 or sub['Mouse_ID'].nunique() < 2:
            return super()._run_single_two_way_anova(column, group_data, design)

        # Fit factorial LMM via R
        result = r_bridge.fit_lmm_factorial(
            sub, 'y ~ factor1 * factor2', '(1|Mouse_ID)', 'factor1', 'factor2',
        )
        if result is None:
            logger.info(f"LMM failed for '{column}', falling back to classical two-way ANOVA")
            return super()._run_single_two_way_anova(column, group_data, design)

        if not result.omnibus:
            logger.warning(f"LMM returned no omnibus for '{column}', falling back to classical")
            return super()._run_single_two_way_anova(column, group_data, design)

        measurement_type = categorize_measurement(column)

        # Store result and data for post-hoc use
        self._lmm_results[column] = result
        self._lmm_sub_data[column] = sub

        # Map omnibus effect names from R to design labels
        # R names: "factor1", "factor2", "factor1:factor2"
        effect_map = {}
        for eff in result.omnibus:
            if eff.effect == 'factor1':
                effect_map['factor1'] = eff.p_value
            elif eff.effect == 'factor2':
                effect_map['factor2'] = eff.p_value
            elif 'factor1' in eff.effect and 'factor2' in eff.effect:
                effect_map['interaction'] = eff.p_value

        # Build reverse mapping: (factor1_value, factor2_value) -> group_name
        factor_to_group = {}
        for g in design.groups:
            f1v = design.factor_mapping[g.name]['factor1']
            f2v = design.factor_mapping[g.name]['factor2']
            factor_to_group[(f1v, f2v)] = g.name

        # Pre-compute simple effects from R emmeans
        simple_effects = []
        for se in result.simple_effects:
            ctx = se.by_level      # level of the conditioning factor
            la, lb = se.level_a, se.level_b  # levels being compared

            # by_factor tells us which factor is held constant
            # so la/lb are levels of the OTHER factor
            if se.by_factor == 'factor2':
                # comparing factor1 levels, factor2 held at ctx
                g1 = factor_to_group.get((la, ctx), la)
                g2 = factor_to_group.get((lb, ctx), lb)
                d1 = sub.loc[(sub['factor1'] == la) & (sub['factor2'] == ctx), 'y']
                d2 = sub.loc[(sub['factor1'] == lb) & (sub['factor2'] == ctx), 'y']
            else:
                # comparing factor2 levels, factor1 held at ctx
                g1 = factor_to_group.get((ctx, la), la)
                g2 = factor_to_group.get((ctx, lb), lb)
                d1 = sub.loc[(sub['factor1'] == ctx) & (sub['factor2'] == la), 'y']
                d2 = sub.loc[(sub['factor1'] == ctx) & (sub['factor2'] == lb), 'y']

            simple_effects.append(StatisticalResult(
                test_name=f"LMM contrast (simple effect, Within {ctx})",
                measurement=column,
                group1_name=g1,
                group1_mean=d1.mean() if len(d1) > 0 else 0.0,
                group1_stderr=d1.std(ddof=1) / math.sqrt(len(d1)) if len(d1) > 1 else 0.0,
                group1_n=len(d1),
                group2_name=g2,
                group2_mean=d2.mean() if len(d2) > 0 else 0.0,
                group2_stderr=d2.std(ddof=1) / math.sqrt(len(d2)) if len(d2) > 1 else 0.0,
                group2_n=len(d2),
                p_value=se.p_value,
                measurement_type=measurement_type,
            ))
        self._lmm_simple_effects[column] = simple_effects

        # Build the 3 omnibus results matching classical format
        f1_levels = sorted(sub['factor1'].unique().tolist())
        f2_levels = sorted(sub['factor2'].unique().tolist())

        results = []
        for p_key, label in [('factor1', design.factor1_name),
                              ('factor2', design.factor2_name),
                              ('interaction', 'Interaction')]:
            p_val = effect_map.get(p_key, 1.0)
            results.append(StatisticalResult(
                test_name=f"{self.name} - {label}",
                measurement=column,
                group1_name=design.groups[0].name,
                group1_mean=0, group1_stderr=0, group1_n=0,
                group2_name=design.groups[1].name,
                group2_mean=0, group2_stderr=0, group2_n=0,
                p_value=p_val,
                measurement_type=measurement_type,
            ))

        # Cache fit info for marginal posthocs (same interface as classical)
        cell_stats_df = sub.groupby(['factor1', 'factor2'])['y'].agg(['mean', 'count']).reset_index()
        cell_stats = {
            (row['factor1'], row['factor2']): {'mean': row['mean'], 'count': int(row['count'])}
            for _, row in cell_stats_df.iterrows()
        }
        self._measurement_fit_info[column] = {
            'mse': None,
            'df_resid': None,
            'cell_stats': cell_stats,
            'factor1_levels': f1_levels,
            'factor2_levels': f2_levels,
        }

        return results

    # ------------------------------------------------------------------
    # Simple effects: return LMM emmeans results or delegate to super()
    # ------------------------------------------------------------------

    def _run_logical_posthoc(self, measurement: str, group_data: Dict[str, pd.DataFrame],
                             design: ExperimentalDesign) -> List[StatisticalResult]:
        """Return pre-computed LMM simple effects or fall back to classical."""
        if measurement in self._lmm_simple_effects:
            return self._lmm_simple_effects[measurement]
        return super()._run_logical_posthoc(measurement, group_data, design)

    # ------------------------------------------------------------------
    # Marginal comparisons: use emmeans marginal or delegate to super()
    # ------------------------------------------------------------------

    def _run_marginal_posthoc(self, measurement: str, group_data: Dict[str, pd.DataFrame],
                              design: ExperimentalDesign, factor_key: str) -> List[StatisticalResult]:
        """Use cached emmeans marginals from the factorial fit, else classical."""

        lmm_result = self._lmm_results.get(measurement)
        sub = self._lmm_sub_data.get(measurement)
        if lmm_result is None or sub is None:
            return super()._run_marginal_posthoc(measurement, group_data, design, factor_key)

        marginal_pairwise = lmm_result.marginals.get(factor_key, [])
        if not marginal_pairwise:
            return super()._run_marginal_posthoc(measurement, group_data, design, factor_key)

        factor_label = design.factor1_name if factor_key == 'factor1' else design.factor2_name
        target_col = factor_key  # 'factor1' or 'factor2'

        measurement_type = categorize_measurement(measurement)

        # Build lookup for model-based emmeans (mean + SE) per level
        emm_lookup = {}
        for emm in lmm_result.emmeans:
            if emm.factor == factor_key:
                emm_lookup[emm.level] = emm

        # Total N per level from raw data
        n_per_level = sub.groupby(target_col)['y'].count().to_dict()

        results = []
        for pw in marginal_pairwise:
            emm_a = emm_lookup.get(pw.level_a)
            emm_b = emm_lookup.get(pw.level_b)

            results.append(StatisticalResult(
                test_name=f"LMM marginal contrast ({factor_label})",
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
