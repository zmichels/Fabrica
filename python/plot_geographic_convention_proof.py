"""Plot a field fabric in eight geographic lower-hemisphere conventions.

The default field case is an included synthetic down-dip example.
Each panel maps a different, explicitly labelled specimen line/pole basis into
that same field fabric. MTEX Python supplies the vectors and rotations; this
script projects their geographic components onto an equal-area stereonet,
matching MTEX's default spherical projection.

Run with the MTEX Python environment, for example::

    python plot_geographic_convention_proof.py

The figure and a machine-readable check report are written in ``python/``.
Pass ``--csv``, ``--sample``, and ``--output`` to select a measured field row.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from mtex import plottingConvention, xvector, yvector, zvector

from fabrica import TP2vec, geoFab


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV = Path(__file__).resolve().parent / "field_fabric_synthetic.csv"
DEFAULT_SAMPLE = "SYNTHETIC DOWN-DIP"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "geographic_convention_proof.png"


@dataclass(frozen=True)
class Case:
    out_label: str
    east_label: str
    out: object
    east: object
    line_label: str
    pole_label: str
    specimen_line: object
    specimen_pole: object


# The first source-axis pair uses +X lineation and +Y foliation pole. The
# other seven are illustrative alternative alignments, kept perpendicular.
CASES = (
    Case("+Z", "+X", zvector, xvector, "+X", "+Y", xvector, yvector),
    Case("+Z", "+Y", zvector, yvector, "+Y", "-X", yvector, -xvector),
    Case("+Z", "-X", zvector, -xvector, "-X", "-Y", -xvector, -yvector),
    Case("+Z", "-Y", zvector, -yvector, "-Y", "+X", -yvector, xvector),
    Case("-Z", "+X", -zvector, xvector, "+X", "-Y", xvector, -yvector),
    Case("-Z", "+Y", -zvector, yvector, "+Y", "+X", yvector, xvector),
    Case("-Z", "-X", -zvector, -xvector, "-X", "+Y", -xvector, yvector),
    Case("-Z", "-Y", -zvector, -yvector, "-Y", "-X", -yvector, -xvector),
)


def xyz(direction) -> np.ndarray:
    """Return Cartesian values from an MTEX Python vector."""
    return np.asarray(direction.data, dtype=float)


def field_row(path: Path, sample: str) -> dict:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row["name"].strip() == sample:
                return {
                    "name": sample,
                    "strike": float(row["strike"]),
                    "dip": float(row["dip"]),
                    "dipDir": row["dipDir"].strip(),
                    "trend": float(row["trend"]),
                    "plunge": float(row["plunge"]),
                }
    raise ValueError(f"No field row named {sample!r} in {path}")


def project_lower(v: np.ndarray, pc: plottingConvention) -> np.ndarray:
    """Equal-area lower-hemisphere coordinates, with geographic E right/N up.

    A unit direction with geographic components (east, north, down) maps to
    (east, north)/sqrt(1 + down), scaled to unit horizon. MTEX's native Schmidt
    plot uses the same projection with a horizon radius of sqrt(2).
    """
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    east = np.sum(v * xyz(pc.east), axis=-1)
    north = np.sum(v * xyz(pc.north), axis=-1)
    down = np.sum(v * xyz(pc.intoScreen), axis=-1)
    if np.any(down < -1e-10):
        raise ValueError("A direction is outside the geographic lower hemisphere")
    return np.stack((east / np.sqrt(1 + down), north / np.sqrt(1 + down)), axis=-1)


def angular_error_degrees(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float) / np.linalg.norm(a)
    b = np.asarray(b, dtype=float) / np.linalg.norm(b)
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b)), np.dot(a, b))))


def compute_case(case: Case, field: dict) -> dict:
    pc = plottingConvention(case.out, case.east)
    strike = field["strike"]
    dip = field["dip"]
    trend = field["trend"]
    plunge = field["plunge"]

    # The supplied plunge is a separate input to the plane-constrained solve.
    measured_line = TP2vec(trend, plunge, pc)
    measured_pole = TP2vec(strike - 90, 90 - dip, pc)
    _, rot, strike_v, solved_pole, solved_line = geoFab(
        strike, dip, trend, case.specimen_line, case.specimen_pole, pc
    )
    rotated_line = rot * case.specimen_line
    rotated_pole = rot * case.specimen_pole

    # A great circle of the input plane runs from one strike direction to
    # its opposite and through down-dip at its middle.
    down_dip = TP2vec(strike + 90, dip, pc)
    theta = np.linspace(0, np.pi, 721)
    trace = (
        np.cos(theta)[:, None] * xyz(strike_v)
        + np.sin(theta)[:, None] * xyz(down_dip)
    )

    plane_miss = float(np.degrees(np.arcsin(np.clip(
        abs(np.dot(xyz(measured_pole), xyz(measured_line))), 0, 1
    ))))
    line_rotation_error = angular_error_degrees(xyz(measured_line), xyz(rotated_line))
    pole_rotation_error = angular_error_degrees(xyz(measured_pole), xyz(rotated_pole))
    solved_line_error = angular_error_degrees(xyz(measured_line), xyz(solved_line))
    solved_pole_error = angular_error_degrees(xyz(measured_pole), xyz(solved_pole))
    trace_plane_error = float(np.max(np.abs(trace @ xyz(measured_pole))))

    if max(plane_miss, line_rotation_error, pole_rotation_error,
           solved_line_error, solved_pole_error) > 1e-8:
        raise ValueError(
            f"{field['name']}: input lineation and foliation or specimen rotation "
            f"are inconsistent for {case.east_label} east / {case.out_label} out "
            f"(maximum angular error {max(plane_miss, line_rotation_error, pole_rotation_error, solved_line_error, solved_pole_error):.6g}°)"
        )
    if trace_plane_error > 1e-10:
        raise AssertionError("Constructed plane trace does not lie in the foliation plane")

    return {
        "convention": f"{case.out_label} out, {case.east_label} east",
        "specimen_line": case.line_label,
        "specimen_pole": case.pole_label,
        "pole_xy": project_lower(xyz(rotated_pole), pc),
        "line_xy": project_lower(xyz(rotated_line), pc),
        "trace_xy": project_lower(trace, pc),
        "plane_miss_degrees": plane_miss,
        "line_rotation_error_degrees": line_rotation_error,
        "pole_rotation_error_degrees": pole_rotation_error,
        "solved_line_error_degrees": solved_line_error,
        "solved_pole_error_degrees": solved_pole_error,
        "trace_plane_max_dot": trace_plane_error,
    }


def check_common_geography(results: list[dict]) -> float:
    first = results[0]
    max_difference = max(
        float(np.max(np.abs(result[key] - first[key])))
        for result in results
        for key in ("pole_xy", "line_xy", "trace_xy")
    )
    if max_difference > 1e-10:
        raise AssertionError(f"Projected geographic plots disagree by {max_difference:g}")
    return max_difference


def draw_stereonet(ax, result: dict) -> None:
    ax.set_aspect("equal")
    ax.set_xlim(-1.17, 1.17)
    ax.set_ylim(-1.17, 1.17)
    ax.axis("off")
    ax.add_patch(Circle((0, 0), 1, facecolor="#fbfdff", edgecolor="#354454", lw=1.2, zorder=0))
    for plunge in (30, 60):
        radius = np.sqrt(1 - np.sin(np.deg2rad(plunge)))
        ax.add_patch(Circle((0, 0), radius, fill=False, edgecolor="#d9e3eb", lw=0.7, ls=":", zorder=1))
    ax.plot((-1, 1), (0, 0), color="#d9e3eb", lw=0.7, ls=":", zorder=1)
    ax.plot((0, 0), (-1, 1), color="#d9e3eb", lw=0.7, ls=":", zorder=1)
    ax.plot(*result["trace_xy"].T, color="#326daa", lw=2.0, zorder=3)
    ax.scatter(*result["pole_xy"], marker="s", s=105, facecolor="#eeb23b",
               edgecolor="#583e10", lw=1.1, zorder=5)
    ax.scatter(*result["line_xy"], marker="o", s=105, facecolor="#d9514e",
               edgecolor="#63201e", lw=1.1, zorder=6)
    for label, x, y in (("N", 0, 1.095), ("E", 1.095, 0),
                        ("S", 0, -1.095), ("W", -1.095, 0)):
        ax.text(x, y, label, ha="center", va="center", fontsize=9.3,
                fontweight="bold", color="#2f3d49")


def render(field: dict, results: list[dict], max_difference: float, output: Path) -> Path:
    fig = plt.figure(figsize=(15.2, 17.2), facecolor="white")
    grid = fig.add_gridspec(
        4, 4, left=0.055, right=0.97, top=0.86, bottom=0.08,
        width_ratios=(1.0, 0.62, 1.0, 0.62), hspace=0.34, wspace=0.16
    )
    fig.suptitle("Same field fabric in eight geographic lower-hemisphere plots",
                 x=0.055, y=0.965, ha="left", fontsize=18, fontweight="bold", color="#243747")
    fig.text(0.055, 0.93,
             f"{field['name']}  •  field geometry: "
             f"strike {field['strike']:g}°, dip {field['dip']:g}° "
             f"({field['dipDir']}), lineation trend {field['trend']:g}°, plunge {field['plunge']:g}°",
             ha="left", fontsize=11.3, color="#344b5d")
    fig.text(0.055, 0.907,
             "Blue = foliation plane trace     ■ = foliation pole     ● = plunging lineation     "
             "N↑  E→  lower hemisphere in every panel",
             ha="left", fontsize=10.5, color="#344b5d")

    for index, result in enumerate(results):
        row, col = divmod(index, 2)
        ax = fig.add_subplot(grid[row, col * 2])
        draw_stereonet(ax, result)
        info = fig.add_subplot(grid[row, col * 2 + 1])
        info.axis("off")
        info.text(0, 0.99, f"{index+1:02d}  {result['convention']}",
                  va="top", fontsize=10, fontweight="bold", color="#243747", transform=info.transAxes)
        info.text(0, 0.82,
                  f"FIELD • {field['name']}\n"
                  f"Strike  {field['strike']:g}°\n"
                  f"Dip     {field['dip']:g}° {field['dipDir']}\n"
                  f"Trend   {field['trend']:g}°\n"
                  f"Plunge  {field['plunge']:g}°",
                  va="top", fontsize=9.4, linespacing=1.4, color="#344b5d", transform=info.transAxes)
        info.text(0, 0.44,
                  "SPECIMEN AXES\n"
                  f"Line    {result['specimen_line']}\n"
                  f"Pole   {result['specimen_pole']}",
                  va="top", fontsize=9.1, linespacing=1.4, color="#344b5d", transform=info.transAxes)
        info.text(0, 0.18,
                  f"On-plane: {result['plane_miss_degrees']:.1e}°\n"
                  f"Rotation: {max(result['line_rotation_error_degrees'], result['pole_rotation_error_degrees']):.1e}°",
                  va="top", fontsize=8.3, linespacing=1.4, color="#486174", transform=info.transAxes)

    fig.text(0.055, 0.047,
             "All eight specimen-axis pairs are illustrative alignments. "
             "Vectors and rotations use MTEX Python; projection is equal-area (Schmidt).",
             ha="left", fontsize=9, color="#566c7b")
    fig.text(0.055, 0.027,
             f"Input lineation lies on the strike/dip plane; "
             f"maximum screen-coordinate difference across panels = {max_difference:.1e}.",
             ha="left", fontsize=9, color="#566c7b")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, facecolor="white")
    pdf_output = output.with_suffix(".pdf")
    fig.savefig(pdf_output, facecolor="white")
    plt.close(fig)
    return pdf_output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="Field case CSV")
    parser.add_argument("--sample", default=DEFAULT_SAMPLE, help="Exact name in the CSV")
    parser.add_argument("--output", type=Path, help="Output PNG path")
    args = parser.parse_args()
    if args.csv is not None and args.output is None:
        parser.error("--output is required with --csv to preserve the bundled synthetic proof")
    args.csv = args.csv or DEFAULT_CSV
    args.output = args.output or DEFAULT_OUTPUT
    if args.output.suffix.lower() != ".png":
        parser.error("--output must end in .png; a matching .pdf and .json are also written")

    field = field_row(args.csv, args.sample)
    results = [compute_case(case, field) for case in CASES]
    max_difference = check_common_geography(results)
    pdf_output = render(field, results, max_difference, args.output)
    report = {
        "source_csv": str(args.csv.relative_to(ROOT)) if args.csv.is_relative_to(ROOT) else str(args.csv),
        "field": field,
        "projection": "equal-area Schmidt, lower hemisphere, geographic N up/E right",
        "max_projected_coordinate_difference": max_difference,
        "cases": [
            {
                key: value.tolist() if isinstance(value, np.ndarray) else value
                for key, value in result.items() if key != "trace_xy"
            }
            for result in results
        ],
    }
    json_output = args.output.with_suffix(".json")
    json_output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PNG: {args.output}\nPDF: {pdf_output}\nChecks: {json_output}")
    print(f"Maximum projected coordinate difference: {max_difference:.3e}")
    print(f"Maximum lineation off-plane angle: {max(row['plane_miss_degrees'] for row in results):.3e}°")
    print(f"Maximum specimen rotation angular error: {max(max(row['line_rotation_error_degrees'], row['pole_rotation_error_degrees']) for row in results):.3e}°")


if __name__ == "__main__":
    main()
