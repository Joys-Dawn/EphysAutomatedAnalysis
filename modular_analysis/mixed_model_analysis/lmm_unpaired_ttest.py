"""
LMM-enhanced unpaired t-test for independent two-group comparisons.

Subclasses UnpairedTTest: overrides the per-measurement test to fit an LMM
with Mouse_ID as a random intercept via R's lmerTest (Satterthwaite df).
Falls back to the classical test (super()) when R is unavailable or the
model fails/is singular.
"""

import os
import math
import logging
import pandas as pd
from typing import Optional

from ..statistical_analysis.tests.unpaired_ttest import UnpairedTTest
from ..shared.data_models import StatisticalResult
from ..shared.utils import categorize_measurement
from .mouse_log import load_mouse_log
from . import r_bridge

logger = logging.getLogger(__name__)


class LMMUnpairedTTest(UnpairedTTest):
    """Unpaired t-test with LMM enhancement for mouse clustering."""

    def __init__(self, mouse_log_path: str):
        super().__init__()
        self.mouse_log = load_mouse_log(mouse_log_path)

    # ------------------------------------------------------------------
    # Data enrichment: add Mouse_ID after normal CSV loading
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
    # Per-measurement override: try LMM via R, fall back to super()
    # ------------------------------------------------------------------

    def _run_single_test(self, column: str, df1: pd.DataFrame, df2: pd.DataFrame,
                         group1_name: str, group2_name: str) -> Optional[StatisticalResult]:
        """Try LMM for one measurement; fall back to classical on failure."""

        if not r_bridge.is_r_available():
            return super()._run_single_test(column, df1, df2, group1_name, group2_name)

        if 'Mouse_ID' not in df1.columns or 'Mouse_ID' not in df2.columns:
            return super()._run_single_test(column, df1, df2, group1_name, group2_name)

        # Build merged long-format frame
        sub1 = df1[['Mouse_ID', column]].dropna(subset=[column]).copy()
        sub1['Group'] = group1_name
        sub2 = df2[['Mouse_ID', column]].dropna(subset=[column]).copy()
        sub2['Group'] = group2_name
        merged = pd.concat([sub1, sub2], ignore_index=True)
        merged.rename(columns={column: 'y'}, inplace=True)

        if len(merged) < 4 or merged['Mouse_ID'].nunique() < 2:
            return super()._run_single_test(column, df1, df2, group1_name, group2_name)

        # Fit LMM via R
        result = r_bridge.fit_lmm_anova(merged, 'y ~ Group', '(1|Mouse_ID)')
        if result is None:
            logger.info(f"LMM failed for '{column}', falling back to classical")
            return super()._run_single_test(column, df1, df2, group1_name, group2_name)

        # Extract p-value (single fixed effect: Group)
        if not result.omnibus:
            logger.warning(f"LMM returned no omnibus for '{column}', falling back to classical")
            return super()._run_single_test(column, df1, df2, group1_name, group2_name)
        p_value = result.omnibus[0].p_value

        # Descriptive stats from raw data (same as classical)
        data1 = df1[column].dropna()
        data2 = df2[column].dropna()

        return StatisticalResult(
            test_name="Linear Mixed Model",
            measurement=column,
            group1_name=group1_name,
            group1_mean=data1.mean(),
            group1_stderr=data1.std(ddof=1) / math.sqrt(len(data1)) if len(data1) > 1 else 0.0,
            group1_n=len(data1),
            group2_name=group2_name,
            group2_mean=data2.mean(),
            group2_stderr=data2.std(ddof=1) / math.sqrt(len(data2)) if len(data2) > 1 else 0.0,
            group2_n=len(data2),
            p_value=p_value,
            measurement_type=categorize_measurement(column),
        )
