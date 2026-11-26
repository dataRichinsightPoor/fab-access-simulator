import json
import math
from pathlib import Path
import unittest
from xml.etree import ElementTree

from fab_access import Binding, Exposure, Profile, compare, simulate, step, summarize, svg_chart, sweep


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((Path(__file__).parents[1] / "example.json").read_text())

    def test_constant_exposure_matches_closed_form(self):
        exposure = Exposure(8, 10)
        for b in (Binding(0.2, 0.4, 0.7), Binding(0.2, 0.4, 0.7, 0.2)):
            r = b.kon_per_nM_h * 8 + b.koff_per_h + b.bound_loss_per_h
            expected = b.accessible_fraction * b.kon_per_nM_h * 8 / r * (1 - math.exp(-r * 4))
            self.assertAlmostEqual(simulate(exposure, b, 4, 0.07, 8)[-1][2], expected, places=12)

    def test_washout_exact_with_no_rebinding(self):
        self.assertAlmostEqual(step(0.6, 0, Binding(1, 0.5, 1, 0.2), 3), 0.6 * math.exp(-2.1))

    def test_accessibility_is_hard_ceiling(self):
        rows = simulate(Exposure(1e8, 20), Binding(1, 0.1, 0.3), 12, 0.05)
        self.assertTrue(all(0 <= r[2] <= 0.3 for r in rows))
        self.assertEqual(summarize(rows, 0.5)["hours_at_or_above_threshold"], 0)

    def test_zero_exposure_gives_zero_occupancy(self):
        self.assertTrue(all(r[2] == 0 for r in simulate(Exposure(0, 1), Binding(0.5, 0.2), 2, 0.1)))

    def test_counterfactual_pairs_share_input(self):
        rows = compare(self.config)
        projection = lambda k: [(t, c) for t, c, _ in rows[k]]
        self.assertEqual(projection("reference"), projection("variant_same_exposure"))
        self.assertEqual(projection("exposure_only"), projection("variant_own_exposure"))
        self.assertEqual(len(list(sweep(self.config))), 42)

    def test_matched_constant_clamp(self):
        self.config["constant_concentration_nM"] = 3
        rows = compare(self.config)
        self.assertEqual(rows["reference"], rows["exposure_only"])
        self.assertEqual(rows["variant_same_exposure"], rows["variant_own_exposure"])

    def test_threshold_crossing_interpolated(self):
        rows = [(0, 1, 0), (2, 1, 1), (4, 1, 0)]
        self.assertEqual(summarize(rows, 0.5)["hours_at_or_above_threshold"], 2)

    def test_midpoint_convergence(self):
        exposure, b = Exposure(10, 5), Binding(0.2, 0.4, 0.8)
        fine = simulate(exposure, b, 12, 0.0125)[-1][2]
        coarse = simulate(exposure, b, 12, 0.2)[-1][2]
        medium = simulate(exposure, b, 12, 0.1)[-1][2]
        self.assertLess(abs(medium - fine), abs(coarse - fine) / 3)

    def test_input_bounds_and_svg(self):
        with self.assertRaises(ValueError):
            Binding(1, 2, 1.1)
        with self.assertRaises(ValueError):
            Exposure(1, 0)
        ElementTree.fromstring(svg_chart(compare(self.config), 0.5, "A & B"))

    def test_custom_profiles_include_both_sets_of_knots_on_a_shared_grid(self):
        self.config["reference_exposure"] = {"points": [[0, 3], [1.007, 3], [72, 0]]}
        self.config["variant_exposure"] = {"points": [[0, 3], [1.013, 1], [72, 4]]}
        output = compare(self.config)
        grid = [t for t, _, _ in output["reference"]]
        self.assertIn(1.007, grid)
        self.assertIn(1.013, grid)
        for rows in output.values():
            self.assertEqual([t for t, _, _ in rows], grid)
        self.assertAlmostEqual(Profile(((0, 3), (2, 1))).at(1), 2)
        with self.assertRaises(ValueError):
            simulate(Profile(((0, 1), (1, 2))), Binding(1, 1), 2, 0.1)


if __name__ == "__main__":
    unittest.main()
