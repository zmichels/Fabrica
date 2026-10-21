"""Deterministic geometric matrix test for Fabrica and MTEX Python.

Run the regression test with the installed MTEX environment::

    python -m unittest -v test_fabrica_matrix

Run this file directly to save a JSON report beside Fabrica::

    python test_fabrica_matrix.py

The oracle below uses independent east/north/down vector formulae. MTEX is
used only for the functions under test and its native equal-area projection.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
import time
import unittest
from pathlib import Path

import numpy as np
from mtex import eareaProjection, plottingConvention, vector3d

from fabrica import V2TP, geoFab


FIBONACCI_COUNT = 48
POLE_ANGLES_DEG = (0, 60, 120, 180, 240, 300)


@dataclass(frozen=True)
class FieldCase:
    name: str
    strike: float
    dip: float
    rake: float


STRIKES_DEG = (0, 89.999, 90, 135, 225, 359.999)
DIPS_DEG = (0, 0.01, 45, 89.99, 90)
RAKES_DEG = (0, 0.01, 45, 90, 135, 179.99, 180)


def _field_grid() -> tuple[tuple[FieldCase, ...], tuple[dict, ...]]:
    """Cross field angles, excluding cases trend alone cannot describe."""
    cases = []
    skipped = []
    for strike in STRIKES_DEG:
        for dip in DIPS_DEG:
            for rake in RAKES_DEG:
                case = FieldCase(f"s{strike:g}-d{dip:g}-r{rake:g}", strike, dip, rake)
                if dip == 90 and rake not in (0, 180):
                    # A vertical foliation's oblique rake is underdetermined by
                    # strike/dip/trend. At rake 90 the trend is undefined.
                    reason = ("undefined_vertical_line_trend" if rake == 90
                              else "vertical_plane_pitch_underdetermined_by_trend")
                    skipped.append({"field": case.name, "reason": reason})
                else:
                    cases.append(case)
    return tuple(cases), tuple(skipped)


FIELD_CASES, SKIPPED_FIELD_CASES = _field_grid()


@dataclass(frozen=True)
class ConventionCase:
    name: str
    east: tuple[float, float, float]
    north: tuple[float, float, float]
    down: tuple[float, float, float]


# These geographic-to-MTEX bases are explicit test data, independent of the
# plottingConvention properties and the TP2vec implementation.
CONVENTIONS = (
    ConventionCase("x-east-z-out", (1, 0, 0), (0, 1, 0), (0, 0, -1)),
    ConventionCase("x-west-z-out", (-1, 0, 0), (0, -1, 0), (0, 0, -1)),
    ConventionCase("x-north-z-out", (0, -1, 0), (1, 0, 0), (0, 0, -1)),
    ConventionCase("x-south-z-out", (0, 1, 0), (-1, 0, 0), (0, 0, -1)),
    ConventionCase("x-east-z-in", (1, 0, 0), (0, -1, 0), (0, 0, 1)),
    ConventionCase("x-west-z-in", (-1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ConventionCase("x-north-z-in", (0, 1, 0), (1, 0, 0), (0, 0, 1)),
    ConventionCase("x-south-z-in", (0, -1, 0), (-1, 0, 0), (0, 0, 1)),
)


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _specimen_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Near-uniform Fibonacci lines plus exact axes; six poles around each."""
    k = np.arange(FIBONACCI_COUNT, dtype=float)
    z = 1 - 2 * (k + 0.5) / FIBONACCI_COUNT
    phi = k * math.pi * (3 - math.sqrt(5))
    radial = np.sqrt(1 - z * z)
    lines = np.column_stack((radial * np.cos(phi), radial * np.sin(phi), z))
    lines = np.vstack((lines, np.eye(3), -np.eye(3)))
    reference = np.where(np.abs(lines[:, 2, None]) < 0.8, (0, 0, 1), (1, 0, 0))
    tangent_a = _unit(np.cross(reference, lines))
    tangent_b = np.cross(lines, tangent_a)
    alpha = np.deg2rad(POLE_ANGLES_DEG)
    poles = (
        np.cos(alpha)[None, :, None] * tangent_a[:, None, :]
        + np.sin(alpha)[None, :, None] * tangent_b[:, None, :]
    )
    line_index = np.repeat(np.arange(len(lines)), len(alpha))
    pole_index = np.tile(np.arange(len(alpha)), len(lines))
    return np.repeat(lines, len(alpha), axis=0), poles.reshape(-1, 3), line_index, pole_index


def _field_oracle(case: FieldCase) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Right-hand strike, down pole, and rake line in geographic ENDown."""
    strike = np.deg2rad(case.strike)
    dip = np.deg2rad(case.dip)
    rake = np.deg2rad(case.rake)
    strike_v = np.array((np.sin(strike), np.cos(strike), 0.0))
    down_dip_v = np.array((np.cos(dip) * np.cos(strike),
                           -np.cos(dip) * np.sin(strike), np.sin(dip)))
    pole_v = np.array((-np.sin(dip) * np.cos(strike),
                        np.sin(dip) * np.sin(strike), np.cos(dip)))
    line_v = np.cos(rake) * strike_v + np.sin(rake) * down_dip_v
    horizontal = np.hypot(line_v[0], line_v[1])
    if horizontal < 1e-8:
        raise AssertionError(f"{case.name} has undefined lineation trend")
    if max(abs(np.linalg.norm(v) - 1) for v in (strike_v, pole_v, line_v)) > 1e-12:
        raise AssertionError(f"nonunit analytic field case: {case.name}")
    if abs(np.dot(line_v, pole_v)) > 1e-12:
        raise AssertionError(f"nonplanar analytic field case: {case.name}")
    trend = float(np.mod(np.rad2deg(np.arctan2(line_v[0], line_v[1])), 360))
    return strike_v, pole_v, line_v, trend


def _equal_area_ens(v: np.ndarray) -> np.ndarray:
    """Independent Schmidt projection, native MTEX horizon radius sqrt(2)."""
    return v[..., :2] * np.sqrt(2 / (1 + v[..., 2]))[..., None]


def _trend_plunge_ens(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    horizontal = np.hypot(v[..., 0], v[..., 1])
    trend = np.mod(np.rad2deg(np.arctan2(v[..., 0], v[..., 1])), 360)
    trend = np.where(horizontal < 1e-10, np.nan, trend)
    plunge = np.rad2deg(np.arctan2(v[..., 2], horizontal))
    plunge = np.where(horizontal < 1e-10, 90, plunge)
    return trend, plunge


def _angular_delta_deg(actual: np.ndarray, expected: np.ndarray) -> np.ndarray:
    both_vertical = np.isnan(actual) & np.isnan(expected)
    mismatch = np.isnan(actual) != np.isnan(expected)
    error = np.abs((actual - expected + 180) % 360 - 180)
    return np.where(mismatch, np.inf, np.where(both_vertical, 0, error))


def _wrong_convention_control() -> dict:
    """Show the direction metric rejects an accidentally swapped plot frame."""
    case = FieldCase("control-southeast-oblique", 135, 45, 45)
    _, _, line_ens, trend = _field_oracle(case)
    expected_xyz = line_ens * np.array((1, 1, -1))
    errors = []
    for convention in CONVENTIONS[:2]:
        east = np.asarray(convention.east, dtype=float)
        down = np.asarray(convention.down, dtype=float)
        pc = plottingConvention(vector3d(-down), vector3d(east))
        _, _, _, _, line = geoFab(
            case.strike, case.dip, trend,
            vector3d((1, 0, 0)), vector3d((0, -1, 0)), pc,
        )
        errors.append(float(np.linalg.norm(np.asarray(line.data) - expected_xyz)))
    return {
        "case": case.name,
        "expected_convention": CONVENTIONS[0].name,
        "incorrect_convention": CONVENTIONS[1].name,
        "correct_direction_error_xyz": errors[0],
        "incorrect_direction_error_xyz": errors[1],
        "detection_threshold_xyz": 1e-8,
        "detected": errors[0] <= 1e-8 and errors[1] > 1e-8,
    }


def _residual_distribution(errors: np.ndarray, *, definition: str = "") -> dict:
    """Summarize raw residuals while retaining a fixed log-scale histogram."""
    values = np.asarray(errors, dtype=float).ravel()
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise AssertionError("residual distribution requires finite nonnegative errors")
    edges = np.linspace(-17.0, -8.0, 19)
    log_values = np.clip(np.log10(np.maximum(values, 1e-17)), edges[0], edges[-1])
    counts, _ = np.histogram(log_values, bins=edges)
    quantiles = np.quantile(values, (0.5, 0.9, 0.99, 0.999))
    return {
        "definition": definition,
        "count": int(values.size),
        "quantiles": {
            "p50": float(quantiles[0]),
            "p90": float(quantiles[1]),
            "p99": float(quantiles[2]),
            "p99_9": float(quantiles[3]),
            "max": float(np.max(values)),
        },
        "log10BinEdges": edges.tolist(),
        "log10BinCounts": counts.tolist(),
        "histogramClipping": "Residuals below 1e-17 enter the first bin; above 1e-8 enter the last bin. Quantiles use raw values.",
    }


def run_matrix() -> dict:
    started = time.perf_counter()
    specimen_lines, specimen_poles, specimen_line_index, pole_angle_index = _specimen_grid()
    n_pairs = len(specimen_lines)
    n_fields = len(FIELD_CASES)
    n_per_convention = n_pairs * n_fields
    source_lines = np.repeat(specimen_lines, n_fields, axis=0)
    source_poles = np.repeat(specimen_poles, n_fields, axis=0)
    source_third = np.cross(source_lines, source_poles)
    field_index = np.tile(np.arange(n_fields), n_pairs)
    pair_index = np.repeat(np.arange(n_pairs), n_fields)
    field_data = [_field_oracle(case) for case in FIELD_CASES]
    field_strikes_ens = np.stack([data[0] for data in field_data])[field_index]
    field_poles_ens = np.stack([data[1] for data in field_data])[field_index]
    field_lines_ens = np.stack([data[2] for data in field_data])[field_index]
    field_trends = np.array([data[3] for data in field_data])[field_index]
    strikes = np.array([case.strike for case in FIELD_CASES])[field_index]
    dips = np.array([case.dip for case in FIELD_CASES])[field_index]

    maxima: dict[str, dict] = {}
    first_projection: dict[str, np.ndarray] = {}
    direction_residuals: list[np.ndarray] = []
    projection_residuals: list[np.ndarray] = []

    def record(metric: str, errors: np.ndarray, convention: str, kind: str = "") -> None:
        errors = np.asarray(errors)
        index = int(np.nanargmax(errors))
        value = float(errors.flat[index])
        if metric in maxima and value <= maxima[metric]["max_error"]:
            return
        pair = int(pair_index[index])
        field = FIELD_CASES[int(field_index[index])]
        maxima[metric] = {
            "max_error": value,
            "worst_case": {
                "convention": convention,
                "component": kind,
                "specimen_line_index": int(specimen_line_index[pair]),
                "specimen_line_xyz": np.round(specimen_lines[pair], 7).tolist(),
                "pole_angle_deg": POLE_ANGLES_DEG[int(pole_angle_index[pair])],
                "field": field.name,
                "strike_deg": field.strike,
                "dip_deg": field.dip,
                "rake_deg": field.rake,
            },
        }

    for convention in CONVENTIONS:
        east = np.asarray(convention.east, dtype=float)
        north = np.asarray(convention.north, dtype=float)
        down = np.asarray(convention.down, dtype=float)
        basis = np.column_stack((east, north, down))
        pc = plottingConvention(vector3d(-down), vector3d(east))
        pc_error = max(
            np.linalg.norm(np.asarray(pc.east.data) - east),
            np.linalg.norm(np.asarray(pc.north.data) - north),
            np.linalg.norm(np.asarray(pc.intoScreen.data) - down),
        )
        record("plotting_basis_xyz", np.full(n_per_convention, pc_error), convention.name)
        expected = {
            "strike": field_strikes_ens @ basis.T,
            "pole": field_poles_ens @ basis.T,
            "lineation": field_lines_ens @ basis.T,
        }
        _, rot, strike_v, pole_v, line_v = geoFab(
            strikes, dips, field_trends,
            vector3d(source_lines), vector3d(source_poles), pc,
        )
        actual = {"strike": np.asarray(strike_v.data),
                  "pole": np.asarray(pole_v.data),
                  "lineation": np.asarray(line_v.data)}
        for kind in ("strike", "pole", "lineation"):
            record("field_direction_xyz", np.linalg.norm(actual[kind] - expected[kind], axis=-1),
                   convention.name, kind)
        mapped = {
            "lineation": np.asarray((rot * vector3d(source_lines)).data),
            "pole": np.asarray((rot * vector3d(source_poles)).data),
            "third_axis": np.asarray((rot * vector3d(source_third)).data),
        }
        mapped_expected = {"lineation": expected["lineation"], "pole": expected["pole"],
                           "third_axis": np.cross(expected["lineation"], expected["pole"])}
        for kind in mapped:
            record("rotation_mapping_xyz", np.linalg.norm(mapped[kind] - mapped_expected[kind], axis=-1),
                   convention.name, kind)
        # Compare both mapped specimen axes in physical geographic ENDown;
        # the basis removes any screen-axis convention before taking norms.
        direction_residuals.append(np.maximum(
            np.linalg.norm(mapped["lineation"] @ basis - field_lines_ens, axis=-1),
            np.linalg.norm(mapped["pole"] @ basis - field_poles_ens, axis=-1),
        ))
        record("pole_plane_dot", np.abs(np.sum(actual["pole"] * actual["lineation"], axis=-1)),
               convention.name)

        projection = eareaProjection(pC=pc, antipodal=False, upper=False)
        mapped_projection_errors = []
        for kind, ens in (("lineation", field_lines_ens), ("pole", field_poles_ens)):
            x, y = projection.project(vector3d(mapped[kind]))
            mapped_xy_unit_horizon = np.stack((x, y), axis=-1) / math.sqrt(2)
            expected_xy_unit_horizon = _equal_area_ens(ens) / math.sqrt(2)
            mapped_projection_errors.append(np.linalg.norm(
                mapped_xy_unit_horizon - expected_xy_unit_horizon, axis=-1,
            ))
        projection_residuals.append(np.maximum(*mapped_projection_errors))
        for kind, v, ens in (("strike", strike_v, field_strikes_ens),
                             ("pole", pole_v, field_poles_ens),
                             ("lineation", line_v, field_lines_ens)):
            x, y = projection.project(v)
            xy = np.stack((x, y), axis=-1)
            oracle_xy = _equal_area_ens(ens)
            record("equal_area_xy", np.linalg.norm(xy - oracle_xy, axis=-1), convention.name, kind)
            if kind in first_projection:
                record("equal_area_cross_convention_xy",
                       np.linalg.norm(xy - first_projection[kind], axis=-1), convention.name, kind)
            else:
                first_projection[kind] = xy

        for kind, v, ens in (("pole", pole_v, field_poles_ens),
                             ("lineation", line_v, field_lines_ens)):
            actual_trend, actual_plunge = V2TP(v, pc)
            expected_trend, expected_plunge = _trend_plunge_ens(ens)
            record("V2TP_trend_deg", _angular_delta_deg(actual_trend, expected_trend),
                   convention.name, kind)
            record("V2TP_plunge_deg", np.abs(actual_plunge - expected_plunge),
                   convention.name, kind)

    thresholds = {
        "plotting_basis_xyz": 1e-10,
        "field_direction_xyz": 1e-8,
        "rotation_mapping_xyz": 1e-8,
        "pole_plane_dot": 1e-8,
        "equal_area_xy": 1e-8,
        "equal_area_cross_convention_xy": 1e-8,
        "V2TP_trend_deg": 1e-7,
        "V2TP_plunge_deg": 1e-7,
    }
    negative_control = _wrong_convention_control()
    return {
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "specimen_line_count": FIBONACCI_COUNT + 6,
        "orthogonal_poles_per_line": len(POLE_ANGLES_DEG),
        "specimen_pair_count": n_pairs,
        "field_grid_strikes_deg": STRIKES_DEG,
        "field_grid_dips_deg": DIPS_DEG,
        "field_grid_rakes_deg": RAKES_DEG,
        "field_grid_candidate_count": len(STRIKES_DEG) * len(DIPS_DEG) * len(RAKES_DEG),
        "field_case_count": n_fields,
        "skipped_field_case_count": len(SKIPPED_FIELD_CASES),
        "skipped_field_cases": SKIPPED_FIELD_CASES,
        "plotting_convention_count": len(CONVENTIONS),
        "total_mapping_cases": n_per_convention * len(CONVENTIONS),
        "thresholds": thresholds,
        "metrics": maxima,
        "residualDistributions": {
            "directionUnitVector": _residual_distribution(
                np.concatenate(direction_residuals),
                definition="Per mapping, maximum Euclidean unit-vector error of rotated specimen line and pole in geographic East/North/Down.",
            ),
            "equalAreaUnitHorizon": _residual_distribution(
                np.concatenate(projection_residuals),
                definition="Per mapping, maximum Euclidean equal-area XY error of rotated specimen line and pole, with horizon radius normalized to one.",
            ),
        },
        "negative_control": negative_control,
        "passed": (all(maxima[k]["max_error"] <= limit for k, limit in thresholds.items())
                   and negative_control["detected"]),
    }


class FabricaMatrixTests(unittest.TestCase):
    def test_residual_distribution_keeps_raw_quantiles_and_counts_every_mapping(self):
        errors = np.array((0.0, 1e-16, 1e-12, 1e-8, 1e-7))
        distribution = _residual_distribution(errors)
        self.assertEqual(distribution["count"], len(errors))
        self.assertEqual(sum(distribution["log10BinCounts"]), len(errors))
        self.assertEqual(distribution["log10BinCounts"][0], 1)
        self.assertEqual(distribution["log10BinCounts"][-1], 2)
        self.assertEqual(distribution["quantiles"]["max"], 1e-7)
        self.assertEqual(distribution["quantiles"]["p50"], 1e-12)
        self.assertEqual(distribution["log10BinEdges"], np.arange(-17, -7.5, 0.5).tolist())

    def test_vertical_line_trend_singularity_is_explicitly_skipped(self):
        singular = [case for case in SKIPPED_FIELD_CASES
                    if case["reason"] == "undefined_vertical_line_trend"]
        ambiguous = [case for case in SKIPPED_FIELD_CASES
                     if case["reason"] == "vertical_plane_pitch_underdetermined_by_trend"]
        self.assertEqual(len(singular), len(STRIKES_DEG))
        self.assertEqual(len(ambiguous), 4 * len(STRIKES_DEG))
        for strike in STRIKES_DEG:
            with self.subTest(strike=strike):
                with self.assertRaisesRegex(AssertionError, "undefined lineation trend"):
                    _field_oracle(FieldCase("vertical-line", strike, 90, 90))

    def test_wrong_plotting_convention_is_detected(self):
        control = _wrong_convention_control()
        self.assertLessEqual(control["correct_direction_error_xyz"],
                             control["detection_threshold_xyz"])
        self.assertGreater(control["incorrect_direction_error_xyz"], 0.5)
        self.assertTrue(control["detected"])

    def test_crossed_specimen_field_convention_matrix(self):
        report = run_matrix()
        self.assertEqual(report["field_grid_candidate_count"], 210)
        self.assertEqual(report["field_case_count"], 180)
        self.assertEqual(report["total_mapping_cases"], 466_560)
        for name in ("directionUnitVector", "equalAreaUnitHorizon"):
            distribution = report["residualDistributions"][name]
            self.assertEqual(distribution["count"], report["total_mapping_cases"])
            self.assertEqual(sum(distribution["log10BinCounts"]), distribution["count"])
            self.assertLessEqual(distribution["quantiles"]["max"], 1e-8)
        for metric, limit in report["thresholds"].items():
            item = report["metrics"][metric]
            self.assertLessEqual(item["max_error"], limit,
                                 f"{metric}: {item['max_error']} at {item['worst_case']}")
        self.assertTrue(report["negative_control"]["detected"])


if __name__ == "__main__":
    report = run_matrix()
    output = Path(__file__).resolve().parent / "geographic_matrix_report_python.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Report: {output}")
    print(f"Mappings: {report['total_mapping_cases']}; passed: {report['passed']}")
    for name, metric in report["metrics"].items():
        print(f"  {name}: {metric['max_error']:.6g}")
    if not report["passed"]:
        raise SystemExit(1)
