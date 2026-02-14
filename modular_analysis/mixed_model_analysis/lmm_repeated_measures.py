"""
LMM-enhanced repeated measures ANOVA for 3+ dependent conditions.

Subclasses RepeatedMeasuresANOVA: fits an LMM with Mouse_ID and Cell_ID
as nested random intercepts via R's lmerTest (Satterthwaite df) for
omnibus, and emmeans for pairwise post-hocs.  Falls back to super()
(classical RM-ANOVA / Friedman) when R is unavailable or the model
fails/is singular.
"""

import os
import math
import logging
import pandas as pd
from typing import Optional, List, Dict

from ..statistical_analysis.tests.repeated_measures_anova import RepeatedMeasuresANOVA
from ..shared.data_models import StatisticalResult
from ..shared.utils import categorize_measurement
from .mouse_log import load_mouse_log
from . import r_bridge

logger = logging.getLogger(__name__)


class LMMRepeatedMeasures(RepeatedMeasuresANOVA):
    """Repeated measures ANOVA with LMM enhancement for mouse clustering."""

    def __init__(self, mouse_log_path: str):
        super().__init__()
        self.name = "Linear Mixed Model (repeated measures)"
        self.mouse_log = load_mouse_log(mouse_log_path)
        # measurement -> list of StatisticalResult (pre-computed pairwise)
        self._lmm_pairwise: Dict[str, List[StatisticalResult]] = {}

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
    # Per-measurement omnibus: try LMM via R, fall back to super()
    # ------------------------------------------------------------------

    def _run_single_rm_anova(self, column: str, all_group_data: Dict[str, pd.DataFrame],
                             conditions: List[str], manifest: pd.DataFrame) -> Optional[StatisticalResult]:
        """Try nested LMM + Satterthwaite; fall back to classical RM-ANOVA."""

        if not r_bridge.is_r_available():
            return super()._run_single_rm_anova(column, all_group_data, conditions, manifest)

        if not all('Mouse_ID' in gdf.columns for gdf in all_group_data.values()):
            return super()._run_single_rm_anova(column, all_group_data, conditions, manifest)

        # Build long-format frame with Cell_ID (= Subject_ID) and Mouse_ID
        long_rows = []
        for condition in conditions:
            if condition not in all_group_data:
                continue
            gdf = all_group_data[condition]
            if column not in gdf.columns:
                continue

            cond_manifest = manifest[manifest['Condition'] == condition]
            fn_to_subject = {
                fn.replace('.abf', ''): sid
                for fn, sid in zip(cond_manifest['Filename'], cond_manifest['Subject_ID'])
            }
            normed = gdf['filename'].str.replace('.abf', '', regex=False)
            sids = normed.map(fn_to_subject)
            mask = sids.notna() & gdf[column].notna() & gdf['Mouse_ID'].notna()

            for sid, val, mid in zip(sids[mask], gdf.loc[mask, column], gdf.loc[mask, 'Mouse_ID']):
                long_rows.append({
                    'Cell_ID': sid,
                    'Condition': condition,
                    'y': val,
                    'Mouse_ID': mid,
                })

        if len(long_rows) < len(conditions) * 3:
            return super()._run_single_rm_anova(column, all_group_data, conditions, manifest)

        sub = pd.DataFrame(long_rows)
        sub['y'] = pd.to_numeric(sub['y'], errors='coerce')
        sub = sub.dropna(subset=['y'])

        # Ensure Cell_ID is unique across mice to avoid crossed RE misspecification
        sub['Cell_ID'] = sub['Mouse_ID'].astype(str) + '_' + sub['Cell_ID'].astype(str)

        # Keep cells present in all conditions
        cell_cond_counts = sub.groupby('Cell_ID')['Condition'].nunique()
        complete_cells = cell_cond_counts[cell_cond_counts == len(conditions)].index
        sub = sub[sub['Cell_ID'].isin(complete_cells)]

        if len(complete_cells) < 3 or sub['Mouse_ID'].nunique() < 2:
            return super()._run_single_rm_anova(column, all_group_data, conditions, manifest)

        # Fit LMM via R with nested RE and emmeans pairwise
        result = r_bridge.fit_lmm_pairwise(
            sub, 'y ~ Condition', '(1|Mouse_ID) + (1|Cell_ID)', 'Condition',
        )
        if result is None:
            logger.info(f"LMM failed for '{column}', falling back to classical RM-ANOVA")
            return super()._run_single_rm_anova(column, all_group_data, conditions, manifest)

        if not result.omnibus:
            logger.warning(f"LMM returned no omnibus for '{column}', falling back to classical")
            return super()._run_single_rm_anova(column, all_group_data, conditions, manifest)
        omnibus_p = result.omnibus[0].p_value

        # Store pairwise for post-hoc phase
        levels = sorted(sub['Condition'].unique().tolist())
        measurement_type = categorize_measurement(column)

        pw_results = []
        for pw in result.pairwise:
            # Descriptive stats from cell-level means
            cell_a = sub.loc[sub['Condition'] == pw.level_a].groupby('Cell_ID')['y'].mean()
            cell_b = sub.loc[sub['Condition'] == pw.level_b].groupby('Cell_ID')['y'].mean()
            common = cell_a.index.intersection(cell_b.index)
            da = cell_a.loc[common]
            db = cell_b.loc[common]
            pw_results.append(StatisticalResult(
                test_name=f"Pairwise LMM ({pw.level_a} vs {pw.level_b})",
                measurement=column,
                group1_name=pw.level_a,
                group1_mean=da.mean() if len(da) > 0 else 0.0,
                group1_stderr=da.std(ddof=1) / math.sqrt(len(da)) if len(da) > 1 else 0.0,
                group1_n=len(da),
                group2_name=pw.level_b,
                group2_mean=db.mean() if len(db) > 0 else 0.0,
                group2_stderr=db.std(ddof=1) / math.sqrt(len(db)) if len(db) > 1 else 0.0,
                group2_n=len(db),
                p_value=pw.p_value,
                measurement_type=measurement_type,
            ))
        self._lmm_pairwise[column] = pw_results

        # Omnibus result (display first two conditions)
        cond_data = {c: sub.loc[sub['Condition'] == c, 'y'] for c in levels}
        c1, c2 = conditions[0], conditions[1]
        d1 = cond_data.get(c1, pd.Series(dtype=float))
        d2 = cond_data.get(c2, pd.Series(dtype=float))

        return StatisticalResult(
            test_name=f"{self.name} (LMM)",
            measurement=column,
            group1_name=c1,
            group1_mean=d1.mean() if len(d1) > 0 else 0.0,
            group1_stderr=d1.std(ddof=1) / math.sqrt(len(d1)) if len(d1) > 1 else 0.0,
            group1_n=len(d1),
            group2_name=c2,
            group2_mean=d2.mean() if len(d2) > 0 else 0.0,
            group2_stderr=d2.std(ddof=1) / math.sqrt(len(d2)) if len(d2) > 1 else 0.0,
            group2_n=len(d2),
            p_value=omnibus_p,
            measurement_type=measurement_type,
        )

    # ------------------------------------------------------------------
    # Pairwise: return stored LMM pairwise or delegate to super()
    # ------------------------------------------------------------------

    def _run_pairwise_comparisons_if_significant(
        self, anova_results: List[StatisticalResult],
        all_group_data: Dict[str, pd.DataFrame],
        conditions: List[str], manifest: pd.DataFrame
    ) -> List[StatisticalResult]:
        """Return LMM emmeans pairwise for LMM measurements, classical for fallbacks."""

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
                classical_anova_results, all_group_data, conditions, manifest
            )

        return lmm_pairwise + classical_pairwise
