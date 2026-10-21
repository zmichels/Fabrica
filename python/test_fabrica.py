"""Behavioral checks against field directions and the MATLAB Fabrica contract."""

from __future__ import annotations

import importlib
import importlib.util
import csv
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np
from mtex import EBSD, calcDensity, crystalFrame, degree, orientation, plottingConvention, rotation, vector3d, xvector, yvector, zvector


def api():
    if importlib.util.find_spec("fabrica") is None:
        raise AssertionError("fabrica.py is required for the Python port")
    return importlib.import_module("fabrica")


def xyz(v):
    return np.asarray(v.data, dtype=float)


CONVENTIONS = (
    # x cardinal, z screen sense, plotting convention, north, east, down
    ("east", "out", plottingConvention(zvector, xvector), (0, 1, 0), (1, 0, 0), (0, 0, -1)),
    ("east", "in", plottingConvention(-zvector, xvector), (0, -1, 0), (1, 0, 0), (0, 0, 1)),
    ("west", "out", plottingConvention(zvector, -xvector), (0, -1, 0), (-1, 0, 0), (0, 0, -1)),
    ("west", "in", plottingConvention(-zvector, -xvector), (0, 1, 0), (-1, 0, 0), (0, 0, 1)),
    ("north", "out", plottingConvention(zvector, -yvector), (1, 0, 0), (0, -1, 0), (0, 0, -1)),
    ("north", "in", plottingConvention(-zvector, yvector), (1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ("south", "out", plottingConvention(zvector, yvector), (-1, 0, 0), (0, 1, 0), (0, 0, -1)),
    ("south", "in", plottingConvention(-zvector, -yvector), (-1, 0, 0), (0, -1, 0), (0, 0, 1)),
)


class DirectionTests(unittest.TestCase):
    def test_identity_cardinals_and_signed_plunge(self):
        f = api()
        np.testing.assert_allclose(xyz(f.TP2vec(0, 0)), (0, 1, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(f.TP2vec(90, 0)), (1, 0, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(f.TP2vec(180, 0)), (0, -1, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(f.TP2vec(270, 0)), (-1, 0, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(f.TP2vec(90, 30)), (np.sqrt(3) / 2, 0, -0.5), atol=1e-12)
        np.testing.assert_allclose(xyz(f.TP2vec(90, -30)), (np.sqrt(3) / 2, 0, 0.5), atol=1e-12)

    def test_rejects_plunge_outside_direction_range(self):
        f = api()
        with self.assertRaisesRegex(ValueError, "Plunge"):
            f.TP2vec(0, 91)
        with self.assertRaisesRegex(ValueError, "Plunge"):
            f.TP2vec(0, -91)

    def test_eight_cardinal_screen_conventions(self):
        f = api()
        for xdir, zdir, convention, north, east, down in CONVENTIONS:
            with self.subTest(xdir=xdir, zdir=zdir):
                np.testing.assert_allclose(xyz(f.TP2vec(0, 0, convention)), north, atol=1e-12)
                np.testing.assert_allclose(xyz(f.TP2vec(90, 0, convention)), east, atol=1e-12)
                np.testing.assert_allclose(xyz(f.TP2vec(0, 90, convention)), down, atol=1e-12)

    def test_inverse_recovers_quadrants_across_conventions(self):
        f = api()
        for xdir, zdir, convention, *_ in CONVENTIONS:
            for trend in (0, 15, 90, 135, 180, 270, 359):
                for plunge in (0, 30, 60):
                    with self.subTest(xdir=xdir, zdir=zdir, trend=trend, plunge=plunge):
                        got_trend, got_plunge = f.V2TP(f.TP2vec(trend, plunge, convention), convention)
                        self.assertAlmostEqual(float(got_trend), trend, places=9)
                        self.assertAlmostEqual(float(got_plunge), plunge, places=9)

    def test_broadcast_array_keeps_shape_and_angles(self):
        f = api()
        trend = np.array([[0, 90], [180, 270]])
        plunge = np.array([10, 50])
        v = f.TP2vec(trend, plunge)
        self.assertEqual(xyz(v).shape, (2, 2, 3))
        got_trend, got_plunge = f.V2TP(v)
        np.testing.assert_allclose(got_trend, trend, atol=1e-9)
        np.testing.assert_allclose(got_plunge, np.broadcast_to(plunge, trend.shape), atol=1e-9)

    def test_vertical_trend_is_undefined(self):
        f = api()
        trend, plunge = f.V2TP(f.TP2vec(37, 90))
        self.assertTrue(np.isnan(trend))
        self.assertAlmostEqual(float(plunge), 90, places=9)

    def test_near_vertical_scaled_vector_uses_relative_tolerance(self):
        f = api()
        trend, plunge = f.V2TP(vector3d(1e-9, 0, -100))
        self.assertTrue(np.isnan(trend))
        self.assertAlmostEqual(float(plunge), 90, places=9)

    def test_upward_vector_is_represented_by_downward_antipode(self):
        f = api()
        trend, plunge = f.V2TP(vector3d(0, 1, 1))
        self.assertAlmostEqual(float(trend), 180, places=9)
        self.assertAlmostEqual(float(plunge), 45, places=9)


class FabricTests(unittest.TestCase):
    def test_trend_rotation_maps_specimen_axes_to_field_fabric(self):
        f = api()
        fab_or, rot, strike, pole, lineation = f.SDT2or(0, 60, 90)
        np.testing.assert_allclose(xyz(strike), (0, 1, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(pole), (-np.sqrt(3) / 2, 0, -0.5), atol=1e-12)
        np.testing.assert_allclose(xyz(lineation), (0.5, 0, -np.sqrt(3) / 2), atol=1e-12)
        np.testing.assert_allclose(xyz(rot * xvector), xyz(lineation), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * -yvector), xyz(pole), atol=1e-9)
        self.assertEqual(type(fab_or).__name__, "orientation")

    def test_trend_strike_parallel_vertical_foliation(self):
        f = api()
        _, rot, strike, pole, lineation = f.SDT2or(0, 90, 0)
        np.testing.assert_allclose(xyz(strike), (0, 1, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(pole), (-1, 0, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(lineation), (0, 1, 0), atol=1e-12)
        np.testing.assert_allclose(xyz(rot * xvector), xyz(lineation), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * -yvector), xyz(pole), atol=1e-9)

    def test_pitch_rotates_strike_line_into_down_dip(self):
        f = api()
        _, _, strike, _, lineation0 = f.SDP2or(0, 60, 0)
        np.testing.assert_allclose(xyz(lineation0), xyz(strike), atol=1e-9)
        _, rot, _, pole, lineation90 = f.SDP2or(0, 60, 90)
        np.testing.assert_allclose(xyz(lineation90), (0.5, 0, -np.sqrt(3) / 2), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * xvector), xyz(lineation90), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * -yvector), xyz(pole), atol=1e-9)

    def test_tiny_upward_pitch_does_not_flip_horizontal_axis(self):
        f = api()
        _, _, _, _, lineation = f.SDP2or(0, 60, -1e-9)
        self.assertGreater(float(lineation.dot(yvector)), 0.999999999)

    def test_into_screen_uses_legacy_specimen_axis_signs(self):
        f = api()
        convention = plottingConvention(-zvector, xvector)
        _, rot, _, pole, lineation = f.SDT2or(0, 60, 90, convention)
        np.testing.assert_allclose(xyz(lineation), (0.5, 0, np.sqrt(3) / 2), atol=1e-9)
        np.testing.assert_allclose(xyz(pole), (-np.sqrt(3) / 2, 0, 0.5), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * -xvector), xyz(lineation), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * yvector), xyz(pole), atol=1e-9)

    def test_explicit_specimen_axes(self):
        f = api()
        _, rot, _, pole, lineation = f.geoFab(0, 60, 90, yvector, -xvector)
        np.testing.assert_allclose(xyz(rot * yvector), xyz(lineation), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * -xvector), xyz(pole), atol=1e-9)

    def test_updip_trend_returns_downward_axial_equivalent(self):
        f = api()
        _, _, _, _, lineation = f.SDT2or(0, 60, 270)
        trend, plunge = f.V2TP(lineation)
        self.assertAlmostEqual(float(trend), 90, places=9)
        self.assertAlmostEqual(float(plunge), 60, places=9)

    def test_rejects_impossible_vertical_plane_trend(self):
        f = api()
        with self.assertRaisesRegex(ValueError, "vertical foliation"):
            f.SDT2or(0, 90, 90)

    def test_rejects_dip_outside_geologic_range(self):
        f = api()
        with self.assertRaisesRegex(ValueError, "Dip"):
            f.SDT2or(0, 100, 90)
        with self.assertRaisesRegex(ValueError, "Dip"):
            f.SDP2or(0, -1, 30)

    def test_rejects_nonperpendicular_specimen_axes(self):
        f = api()
        with self.assertRaisesRegex(ValueError, "perpendicular"):
            f.geoFab(0, 60, 90, xvector, xvector)

    def test_field_fabric_quadrants_across_eight_conventions(self):
        f = api()
        for xdir, zdir, convention, north, east, down in CONVENTIONS:
            with self.subTest(xdir=xdir, zdir=zdir):
                _, _, strike, pole, lineation = f.SDT2or(0, 60, 90, convention)
                _, _, _, _, pitched = f.SDP2or(0, 60, 90, convention)
                np.testing.assert_allclose(xyz(strike), north, atol=1e-9)
                expected_pole = -np.sqrt(3) / 2 * np.asarray(east) + 0.5 * np.asarray(down)
                expected_line = 0.5 * np.asarray(east) + np.sqrt(3) / 2 * np.asarray(down)
                np.testing.assert_allclose(xyz(pole), expected_pole, atol=1e-9)
                np.testing.assert_allclose(xyz(lineation), expected_line, atol=1e-9)
                np.testing.assert_allclose(xyz(pitched), expected_line, atol=1e-9)

    def test_array_fabric_rotations_map_each_case(self):
        f = api()
        strikes = np.array([0, 90])
        dips = np.array([60, 45])
        trends = np.array([90, 180])
        _, rot, _, pole, lineation = f.SDT2or(strikes, dips, trends)
        self.assertEqual(rot.shape, (2,))
        np.testing.assert_allclose(xyz(rot * xvector), xyz(lineation), atol=1e-9)
        np.testing.assert_allclose(xyz(rot * -yvector), xyz(pole), atol=1e-9)

    def test_dip_direction_rake_uses_distinct_specimen_source_axes(self):
        f = api()
        fab_or, strikes, poles, lines = f.FandL2Or(90, 60, 90)
        np.testing.assert_allclose(xyz(strikes), (0, 1, 0), atol=1e-9)
        np.testing.assert_allclose(xyz(poles), (-np.sqrt(3) / 2, 0, -0.5), atol=1e-9)
        np.testing.assert_allclose(xyz(lines), (0.5, 0, -np.sqrt(3) / 2), atol=1e-9)
        np.testing.assert_allclose(xyz(fab_or * zvector), xyz(poles), atol=1e-9)
        np.testing.assert_allclose(xyz(fab_or * xvector), xyz(lines), atol=1e-9)

    def test_dip_direction_rake_accepts_arrays(self):
        f = api()
        dip_direction = np.array([90, 180, 270])
        dip = np.array([30, 45, 60])
        rake = np.array([0, 45, 90])
        fab_or, strikes, poles, lines = f.FandL2Or(dip_direction, dip, rake)
        self.assertEqual(fab_or.shape, (3,))
        np.testing.assert_allclose(xyz(fab_or * zvector), xyz(poles), atol=1e-9)
        np.testing.assert_allclose(xyz(fab_or * xvector), xyz(lines), atol=1e-9)
        np.testing.assert_allclose(np.sum(xyz(lines) * xyz(poles), axis=-1), 0, atol=1e-9)

    def test_synthetic_example_csv_all_field_pairs(self):
        f = api()
        path = Path(__file__).resolve().parent / "example.csv"
        with path.open(newline="") as source:
            rows = list(csv.DictReader(source))
        self.assertGreaterEqual(len(rows), 10)
        dd = np.array([float(row["dipDir"]) for row in rows])
        dip = np.array([float(row["dip"]) for row in rows])
        rake = np.array([float(row["rhrake"]) for row in rows])
        fab_or, _, poles, lines = f.FandL2Or(dd, dip, rake)
        self.assertEqual(fab_or.shape, (len(rows),))
        np.testing.assert_allclose(xyz(fab_or * zvector), xyz(poles), atol=1e-8)
        np.testing.assert_allclose(xyz(fab_or * xvector), xyz(lines), atol=1e-8)
        np.testing.assert_allclose(np.sum(xyz(lines) * xyz(poles), axis=-1), 0, atol=1e-8)


class PlotOutputSafetyTests(unittest.TestCase):
    def test_custom_csv_requires_explicit_output(self):
        here = Path(__file__).resolve().parent
        csv_path = here / "field_fabric_synthetic.csv"
        for script in ("plot_geographic_convention_proof.py", "plot_geographic_matrix_figure.py"):
            with self.subTest(script=script):
                result = subprocess.run(
                    [sys.executable, str(here / script), "--csv", str(csv_path)],
                    text=True, capture_output=True, check=False,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("--output is required with --csv", result.stderr)


class DataRotationTests(unittest.TestCase):
    def setUp(self):
        self.cs = crystalFrame("mmm", mineral="TestF")
        self.ori = orientation(rotation.id, self.cs)
        self.ebsd = EBSD(
            np.array([[0.0, 0.0], [1.0, 0.0]]),
            rotation([[1, 0, 0, 0], [1, 0, 0, 0]]),
            np.array([1, 1]),
            [self.cs],
        )
        self.rot = rotation.byAxisAngle(zvector, 90 * degree)

    def test_orientation_rotation_is_left_composition(self):
        f = api()
        moved = f.rotTex2Geo(self.ori, self.rot)
        np.testing.assert_allclose(xyz(moved * xvector), xyz(self.rot * xvector), atol=1e-9)
        np.testing.assert_allclose(xyz(self.ori * xvector), xyz(xvector), atol=1e-9)

    def test_ebsd_rotation_preserves_map_positions(self):
        f = api()
        moved = f.rotEBSD2Geo(self.ebsd, self.rot)
        np.testing.assert_allclose(moved.pos.data, self.ebsd.pos.data, atol=1e-12)
        np.testing.assert_allclose(xyz(moved.orientations * xvector), np.broadcast_to(xyz(yvector), (2, 3)), atol=1e-9)
        np.testing.assert_allclose(xyz(self.ebsd.orientations * xvector), np.broadcast_to(xyz(xvector), (2, 3)), atol=1e-9)

    def test_odf_rotation_moves_density_peak(self):
        f = api()
        odf = calcDensity(self.ori, halfwidth=10 * degree)
        moved = f.rotODFGeo(odf, self.rot)
        self.assertAlmostEqual(float(moved.eval(self.rot * self.ori)), float(odf.eval(self.ori)), places=6)
        self.assertLess(float(moved.eval(self.ori)), float(moved.eval(self.rot * self.ori)))

    def test_odf_rotation_accepts_sample_dictionary(self):
        f = api()
        sample = f.makeGeoSample(0, 60, 90)
        odf = calcDensity(self.ori, halfwidth=10 * degree)
        moved = f.rotODFGeo(odf, sample)
        expected_peak = sample["geo"]["fabRot"] * self.ori
        self.assertAlmostEqual(float(moved.eval(expected_peak)), float(odf.eval(self.ori)), places=6)

    def test_make_geographic_sample_rotates_map_and_named_phase(self):
        f = api()
        sample = f.makeGeoSample(0, 60, 90, self.ebsd, "TestF")
        self.assertEqual(sample["geo"]["strike"], 0)
        self.assertEqual(sample["geo"]["dip"], 60)
        self.assertEqual(sample["geo"]["trend"], 90)
        self.assertEqual(len(sample["spec"]["ebsd"]), 1)
        np.testing.assert_allclose(sample["geo"]["rotebsd"][0].pos.data, self.ebsd.pos.data, atol=1e-12)
        moved = sample["geo"]["TestF"]
        expected = np.broadcast_to(xyz(sample["geo"]["fabRot"] * xvector), (2, 3))
        np.testing.assert_allclose(xyz(moved * xvector), expected, atol=1e-9)

    def test_sample_rake_can_exceed_ninety_and_plunge_is_unrounded(self):
        f = api()
        sample = f.makeGeoSample(0, 60, 135)
        self.assertAlmostEqual(float(sample["geo"]["rake"]), 116.5650511771, places=7)
        self.assertAlmostEqual(float(sample["geo"]["plunge"]), 50.7684795164, places=7)

    def test_sample_accepts_positional_convention_and_stores_it(self):
        f = api()
        pc = plottingConvention(-zvector, xvector)
        sample = f.makeGeoSample(0, 60, 90, pc)
        self.assertEqual(sample["geo"]["how2plot"], pc)
        self.assertEqual(sample["spec"]["how2plot"], pc)
        np.testing.assert_allclose(xyz(sample["spec"]["lin"]), xyz(-xvector), atol=1e-12)

    def test_sample_records_measured_and_downward_lineation_trends(self):
        f = api()
        sample = f.makeGeoSample(0, 60, 270)
        self.assertEqual(sample["geo"]["trend"], 270)
        self.assertAlmostEqual(float(sample["geo"]["lineTrend"]), 90, places=9)
        self.assertAlmostEqual(float(sample["geo"]["plunge"]), 60, places=9)


if __name__ == "__main__":
    unittest.main()
