"""
LMM-enhanced paired t-test for dependent two-group comparisons.

Subclasses PairedTTest: fits an LMM with Mouse_ID and Cell_ID as nested
random intercepts via R's lmerTest (Satterthwaite df).  Falls back to
super() (classical paired t-test) when R is unavailable or the model
fails/is singular.
"""

import os
import math
import logging
import pandas as pd
from typing import Optional

from ..statistical_analysis.tests.paired_ttest import PairedTTest
from ..shared.data_models import StatisticalResult
from ..shared.utils import categorize_measurement
from .mouse_log import load_mouse_log
from . import r_bridge

logger = logging.getLogger(__name__)


class LMMPairedTTest(PairedTTest):
    """Paired t-test with LMM enhancement for mouse clustering."""

    def __init__(self, mouse_log_path: str):
        super().__init__()
        self.mouse_log = load_mouse_log(mouse_log_path)

    # ------------------------------------------------------------------
    # Data enrichment: add Mouse_ID after parent builds unified_df
    # ------------------------------------------------------------------

    def _create_unified_dataframe(self, df1, df2, manifest, cond1_name, cond2_name):
        """Add Mouse_ID column after the normal unified dataframe creation."""
        unified = super()._create_unified_dataframe(df1, df2, manifest, cond1_name, cond2_name)
        if 'filename' in unified.columns:
            unified['Mouse_ID'] = unified['filename'].apply(
                lambda f: self.mouse_log.get(os.path.basename(str(f).strip()))
            )
            missing = unified['Mouse_ID'].isna().sum()
            if missing > 0:
                logger.warning(f"Paired data: {missing}/{len(unified)} files not in mouse log, dropping")
                unified = unified.dropna(subset=['Mouse_ID'])
        return unified

    # ------------------------------------------------------------------
    # Per-measurement override: try LMM with nested RE, fall back
    # ------------------------------------------------------------------

    def _run_single_paired_test(self, column: str, unified_df: pd.DataFrame,
                                cond1_name: str, cond2_name: str) -> Optional[StatisticalResult]:
        """Try nested LMM (1|Mouse_ID) + (1|Cell_ID); fall back to classical."""

        if not r_bridge.is_r_available():
            return super()._run_single_paired_test(column, unified_df, cond1_name, cond2_name)

        if 'Mouse_ID' not in unified_df.columns:
            return super()._run_single_paired_test(column, unified_df, cond1_name, cond2_name)

        sub = unified_df[['Subject_ID', 'Condition', 'Mouse_ID', column]].dropna(subset=[column]).copy()

        # Build long-format with Cell_ID = Subject_ID
        # Ensure Cell_ID is unique across mice to avoid crossed RE misspecification
        sub = sub.rename(columns={'Subject_ID': 'Cell_ID', column: 'y'})
        sub['Cell_ID'] = sub['Mouse_ID'].astype(str) + '_' + sub['Cell_ID'].astype(str)
        sub['y'] = pd.to_numeric(sub['y'], errors='coerce')
        sub = sub.dropna(subset=['y'])

        # Keep only cells present in both conditions
        cell_cond = sub.groupby('Cell_ID')['Condition'].nunique()
        complete = cell_cond[cell_cond == 2].index
        sub = sub[sub['Cell_ID'].isin(complete)]

        if len(complete) < 2 or sub['Mouse_ID'].nunique() < 2:
            return super()._run_single_paired_test(column, unified_df, cond1_name, cond2_name)

        # Fit LMM via R with nested RE
        result = r_bridge.fit_lmm_anova(
            sub, 'y ~ Condition', '(1|Mouse_ID) + (1|Cell_ID)',
        )
        if result is None:
            logger.info(f"LMM failed for '{column}', falling back to classical paired test")
            return super()._run_single_paired_test(column, unified_df, cond1_name, cond2_name)

        if not result.omnibus:
            logger.warning(f"LMM returned no omnibus for '{column}', falling back to classical")
            return super()._run_single_paired_test(column, unified_df, cond1_name, cond2_name)
        p_value = result.omnibus[0].p_value

        # Descriptive stats from cell-level means
        cell_means = sub.groupby(['Cell_ID', 'Condition'])['y'].mean().reset_index()
        data1 = cell_means.loc[cell_means['Condition'] == cond1_name, 'y']
        data2 = cell_means.loc[cell_means['Condition'] == cond2_name, 'y']

        return StatisticalResult(
            test_name="Linear Mixed Model (paired)",
            measurement=column,
            group1_name=cond1_name,
            group1_mean=data1.mean() if len(data1) > 0 else 0.0,
            group1_stderr=data1.std(ddof=1) / math.sqrt(len(data1)) if len(data1) > 1 else 0.0,
            group1_n=len(data1),
            group2_name=cond2_name,
            group2_mean=data2.mean() if len(data2) > 0 else 0.0,
            group2_stderr=data2.std(ddof=1) / math.sqrt(len(data2)) if len(data2) > 1 else 0.0,
            group2_n=len(data2),
            p_value=p_value,
            measurement_type=categorize_measurement(column),
        )
