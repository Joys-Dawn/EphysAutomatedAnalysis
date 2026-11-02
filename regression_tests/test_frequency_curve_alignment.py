"""Regression tests for whole-curve models using complete recordings.

Run with: python -B -m unittest discover -s regression_tests -v
"""

import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from modular_analysis.data_extraction.extractor import DataExtractor
from modular_analysis.shared.config import AnalysisConfig, ProtocolConfig
from modular_analysis.shared.data_models import DesignType, ExperimentalDesign, GroupInfo
from modular_analysis.statistical_analysis.designs import DesignManager
from modular_analysis.statistical_analysis.tests import frequency_analysis as independent
from modular_analysis.statistical_analysis.tests import frequency_analysis_dependent as dependent


class FrequencyCurveAlignmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "Results").mkdir()
        self.config = ProtocolConfig()
        self.independent = independent.FrequencyAnalyzer(self.config)
        self.dependent = dependent.FrequencyAnalyzerDependent(self.config)

    def make_raw(self, group, offset=0):
        currents = np.arange(-60, 441, 20, dtype=float)
        columns = {}
        for cell, rheobase in enumerate((20, 40, 60)):
            frequencies = np.where(currents >= rheobase, offset + 100 * (cell + 1) + np.arange(len(currents)), 0.0)
            # A zero after firing begins is still an observation.
            frequencies[10] = 0
            columns[f"{group}_{cell}.abf"] = currents
            columns[f"Values_{cell}"] = frequencies
        raw = pd.DataFrame(columns)
        self.assertFalse(raw.isna().any().any())
        path = self.root / "Results" / f"Calc_{group}_frequency_vs_current.csv"
        raw.to_csv(path, index=False)
        return path, raw

    def expected_fold_rows(self, raw, group, integer_only=False):
        rows = []
        currents = self.config.min_current + np.arange(len(raw)) * self.config.step_size
        for cell, column in enumerate(raw.filter(regex="^Values_").columns):
            values = raw[column].to_numpy()
            first = np.flatnonzero(values > 0)[0]
            for index in range(first, len(raw)):
                fold = currents[index] / currents[first]
                if fold > 10:
                    continue
                if integer_only and fold not in range(1, 11):
                    continue
                rows.append({"Cell_id": f"{group}_{cell}", "FoldRheobase": fold, "Frequency": values[index]})
        return rows

    def capture_independent(self, data, analysis, design=None):
        # Stop at the model boundary so assertions inspect real preparation code.
        with patch.object(independent, "mixedlm", side_effect=RuntimeError("capture model input")) as model:
            with self.assertLogs(independent.logger, level="ERROR"):
                self.independent._run_unified_mixed_effects_frequency(data, analysis, design)
        self.assertEqual(model.call_count, 1)
        return model.call_args.args[1]

    def capture_paired(self, data, analysis, design):
        with patch.object(dependent.smf, "mixedlm", side_effect=RuntimeError("capture model input")) as model:
            with self.assertLogs(dependent.logger, level="ERROR"):
                self.dependent._run_global_mixed_effects_paired(data, analysis, ["A", "B"], design)
        self.assertEqual(model.call_count, 1)
        return model.call_args.args[1]

    def test_fold_table_preserves_cells_folds_and_export(self):
        path, raw = self.make_raw("A")
        table = self.independent._make_rheo_df(str(path), "A", str(self.root), integers_only=False)
        table = table.set_index("Fold Rheobase")
        self.assertTrue(table.index.is_monotonic_increasing)
        self.assertEqual(list(table.filter(regex="^Values_").columns), ["Values_0", "Values_1", "Values_2"])
        for expected in self.expected_fold_rows(raw, "A"):
            column = "Values_" + expected["Cell_id"].rsplit("_", 1)[1]
            self.assertEqual(table.loc[expected["FoldRheobase"], column], expected["Frequency"])
        self.assertTrue(pd.isna(table.loc[1.5, "Values_0"]))
        # The same complete current protocol spans fewer folds at 60 pA rheobase.
        self.assertTrue(pd.isna(table.loc[10, "Values_2"]))
        self.assertEqual(table.loc[1, "n"], 3)
        self.assertEqual(table.loc[1, "mean"], (104 + 205 + 306) / 3)
        exported = pd.read_csv(self.root / "Results" / "Calc_A_frequency_vs_fold_rheobase.csv")
        for cell, first in enumerate((4, 5, 6)):
            np.testing.assert_allclose(exported[f"Values_{cell}"], raw[f"Values_{cell}"].shift(-first), equal_nan=True)
            expected_folds = raw[f"A_{cell}.abf"].shift(-first) / raw[f"A_{cell}.abf"].iloc[first]
            np.testing.assert_allclose(exported[f"A_{cell}.abf"], expected_folds, equal_nan=True)

    def test_integer_folds_keep_the_correct_observations_and_summaries(self):
        path, raw = self.make_raw("A")
        table = self.independent._make_rheo_df(str(path), integers_only=True).set_index("Fold Rheobase")
        expected = pd.DataFrame(self.expected_fold_rows(raw, "A", integer_only=True))
        self.assertTrue(set(table.index).issubset(range(1, 11)))
        for fold, rows in expected.groupby("FoldRheobase"):
            self.assertEqual(table.loc[fold, "n"], len(rows))
            self.assertAlmostEqual(table.loc[fold, "mean"], rows.Frequency.mean())
            if len(rows) > 1:
                self.assertAlmostEqual(table.loc[fold, "se"], rows.Frequency.sem())

    def test_whole_curve_fold_inputs_for_independent_and_factorial_designs(self):
        for design_type, group_count in ((DesignType.INDEPENDENT_TWO_GROUP, 2), (DesignType.INDEPENDENT_MULTI_GROUP, 3), (DesignType.FACTORIAL_2X2, 4)):
            with self.subTest(design=design_type):
                data, expected = {}, []
                names = [f"G{i}" for i in range(group_count)]
                for index, name in enumerate(names):
                    path, raw = self.make_raw(name, offset=index * 1000)
                    table = self.independent._make_rheo_df(str(path), integers_only=False)
                    data[name] = table[table["Fold Rheobase"] <= 10]
                    expected.extend(self.expected_fold_rows(raw, name))
                design = ExperimentalDesign(design_type, [GroupInfo(name, "") for name in names], factor_mapping={name: {"factor1": str(i // 2), "factor2": str(i % 2)} for i, name in enumerate(names)})
                actual = self.capture_independent(data, "fold_rheobase", design)
                columns = ["Cell_id", "FoldRheobase", "Frequency"]
                pd.testing.assert_frame_equal(actual[columns].sort_values(columns).reset_index(drop=True), pd.DataFrame(expected)[columns].sort_values(columns).reset_index(drop=True))

    def test_paired_complete_manifest_controls_pairing_despite_filename_order(self):
        manifest = pd.DataFrame({"Subject ID": ["S1", "S2"], "Group": ["G", "G"], "A": ["A_0001.abf", "A_0002.abf"], "B": ["B_0002.abf", "B_0001.abf"]})
        with patch.object(pd, "read_excel", return_value=manifest):
            design = DesignManager.create_paired_two_group(GroupInfo("A", ""), GroupInfo("B", ""), "in_memory_fixture.xlsx")
        lookup = dict(zip(design.pairing_manifest.Filename, design.pairing_manifest.Subject_ID))
        currents = [-60, -40, -20, 0, 20, 40, 60, 80]
        raw = {
            "A": pd.DataFrame({"A_0001.abf": currents, "Values_0": [0, 0, 0, 0, 11, 12, 13, 14], "A_0002.abf": currents, "Values_1": [0, 0, 0, 0, 0, 21, 22, 23]}, dtype=float),
            "B": pd.DataFrame({"B_0001.abf": currents, "Values_0": [0, 0, 0, 0, 0, 31, 32, 33], "B_0002.abf": currents, "Values_1": [0, 0, 0, 0, 41, 42, 43, 44]}, dtype=float),
        }
        extractor = DataExtractor(AnalysisConfig())
        data = {}
        expected = []
        for group, frame in raw.items():
            self.assertFalse(frame.isna().any().any())
            path = self.root / "Results" / f"Calc_{group}_frequency_vs_current.csv"
            frame.to_csv(path, index=False)
            extractor._add_subject_id_to_wide_csv(str(path.parent), group, "frequency_vs_current", lookup)
            data[group] = pd.read_csv(path, index_col=False)
            for cell in range(2):
                filename = f"{group}_{cell + 1:04d}.abf"
                values = frame[f"Values_{cell}"].to_numpy()
                rheobase = currents[np.flatnonzero(values > 0)[0]]
                for current, frequency in zip(currents, values):
                    if current > 0:
                        expected.append({"Subject_ID": "Subject_" + lookup[filename], "Condition": group, "Current": float(current), "FoldRheobase": current / rheobase, "Frequency": float(frequency)})
        for analysis, step in (("current", "Current"), ("fold_rheobase", "FoldRheobase")):
            with self.subTest(analysis=analysis):
                actual = self.capture_paired(data, analysis, design)
                reference = pd.DataFrame(expected)
                if analysis == "fold_rheobase":
                    reference = reference[reference.FoldRheobase >= 1]
                columns = ["Subject_ID", "Condition", step, "Frequency"]
                pd.testing.assert_frame_equal(actual[columns].sort_values(columns).reset_index(drop=True), reference[columns].sort_values(columns).reset_index(drop=True))

    def test_paired_p_values_do_not_depend_on_recording_order(self):
        rng = np.random.default_rng(20261029)
        currents = np.arange(-60, 401, 20)
        folds = currents[4:] / 20
        curves = {}
        subjects = list(range(12))
        for subject in subjects:
            intercept, slope = rng.normal(0, 5), rng.normal(0, 0.18)
            for ci, condition in enumerate(("A", "B")):
                positive = 25 + intercept + (0.9 + slope) * (folds - 10.5) + ci * (0.4 + 0.03 * (folds - 10.5)) + rng.normal(0, 1, len(folds))
                self.assertTrue(np.all(positive > 0))
                curves[condition, subject] = np.r_[np.zeros(4), positive]

        def wide(condition, order):
            columns = {}
            for position, subject in enumerate(order):
                columns[f"Subject_S{subject}"] = currents
                columns[f"Values_{position}"] = curves[condition, subject]
            return pd.DataFrame(columns)

        design = ExperimentalDesign(DesignType.PAIRED_TWO_GROUP, [], within_factor_name="Condition")
        for analysis in ("current", "fold_rheobase"):
            with self.subTest(analysis=analysis):
                results = []
                for order in (subjects, subjects[::-1]):
                    data = {"A": wide("A", subjects), "B": wide("B", order)}
                    with warnings.catch_warnings():
                        warnings.simplefilter("error")
                        result = self.dependent._run_global_mixed_effects_paired(data, analysis, ["A", "B"], design)
                    self.assertIsNotNone(result)
                    results.append([row["p-value"] for row in result])
                np.testing.assert_allclose(results[0], results[1], rtol=1e-8, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
