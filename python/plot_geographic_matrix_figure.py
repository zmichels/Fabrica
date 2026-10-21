"""Illustrate geographic agreement and the full matrix residual distributions.

The stereonets overlay eight MTEX Python plotting/specimen conventions. Their
great circles and symbols are calculated by the existing proof code. The two
density panels read per-mapping residual histograms from the independently run
MATLAB and Python matrix audits. Run the audits first, then this script.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np

plt.rcParams["mathtext.fontset"] = "dejavusans"

from plot_geographic_convention_proof import (
    CASES,
    DEFAULT_CSV,
    DEFAULT_SAMPLE,
    check_common_geography,
    compute_case,
    draw_stereonet,
    field_row,
)
from test_fabrica_matrix import FieldCase, _field_oracle


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "geographic_matrix_figure.png"
MATLAB_REPORT = ROOT / "matlab" / "mtex61" / "geographic_matrix_report_mtex61.json"
PYTHON_REPORT = Path(__file__).resolve().parent / "geographic_matrix_report_python.json"
COLORS = {"MATLAB MTEX 6.1": "#176178", "MTEX Python": "#b35a36"}


def load_reports(matlab_path: Path, python_path: Path) -> dict:
    reports = {
        "MATLAB MTEX 6.1": json.loads(matlab_path.read_text(encoding="utf-8")),
        "MTEX Python": json.loads(python_path.read_text(encoding="utf-8")),
    }
    for name, report in reports.items():
        total = report["evaluations"] if name.startswith("MATLAB") else report["total_mapping_cases"]
        for metric in ("directionUnitVector", "equalAreaUnitHorizon"):
            distribution = report["residualDistributions"][metric]
            if distribution["count"] != total:
                raise ValueError(f"{name}: {metric} count does not match mapping count")
            if sum(distribution["log10BinCounts"]) != total:
                raise ValueError(f"{name}: {metric} histogram does not sum to mapping count")
            if len(distribution["log10BinEdges"]) != len(distribution["log10BinCounts"]) + 1:
                raise ValueError(f"{name}: {metric} histogram edges/counts mismatch")
        if name == "MTEX Python" and not report["passed"]:
            raise ValueError("Python matrix audit did not pass")
    return reports


def oblique_field() -> dict:
    # This 225/45/45 synthetic field case is present in both audit grids.
    case = FieldCase("synthetic-oblique", 225, 45, 45)
    _, _, line, trend = _field_oracle(case)
    return {
        "name": "Synthetic oblique grid case",
        "strike": case.strike,
        "dip": case.dip,
        "dipDir": "NW",
        "trend": trend,
        "plunge": float(np.degrees(np.arcsin(line[2]))),
    }


def plotted_examples(csv_path: Path, sample: str) -> list[tuple[dict, list[dict], float]]:
    fields = (field_row(csv_path, sample), oblique_field())
    output = []
    for field in fields:
        cases = [compute_case(case, field) for case in CASES]
        output.append((field, cases, check_common_geography(cases)))
    return output


def overlay_stereonet(ax, cases: list[dict]) -> None:
    draw_stereonet(ax, cases[0])
    # Draw every additional convention at its own projected coordinates. Their
    # sub-pixel differences deliberately remain invisible at figure scale.
    for case in cases[1:]:
        ax.plot(*case["trace_xy"].T, color="#326daa", lw=1.8, alpha=0.16, zorder=3)
        ax.scatter(*case["pole_xy"], marker="s", s=105, facecolor="#eeb23b",
                   edgecolor="#583e10", lw=1.0, alpha=0.16, zorder=5)
        ax.scatter(*case["line_xy"], marker="o", s=105, facecolor="#d9514e",
                   edgecolor="#63201e", lw=1.0, alpha=0.16, zorder=6)


def density_panel(ax, reports: dict, metric: str, title: str) -> None:
    for name, report in reports.items():
        distribution = report["residualDistributions"][metric]
        edges = np.power(10.0, np.asarray(distribution["log10BinEdges"], dtype=float))
        fractions = np.asarray(distribution["log10BinCounts"], dtype=float) / distribution["count"]
        ax.stairs(fractions, edges, color=COLORS[name], lw=2.1, label=name)
    ax.set_xscale("log")
    ax.set_xlim(1e-17, 1e-8)
    ax.set_xticks([1e-17, 1e-15, 1e-13, 1e-11, 1e-9])
    ax.set_ylabel("Fraction per 0.5 decade", fontsize=9.5)
    ax.set_xlabel("Residual magnitude", fontsize=9.5)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    ax.grid(axis="y", color="#dae3e9", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=11.5, fontweight="bold", color="#243747", loc="left", pad=9)
    ax.tick_params(axis="both", labelsize=8.7)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    maxima = [reports[name]["residualDistributions"][metric]["quantiles"]["max"]
              for name in reports]
    ax.text(0.01, 0.97, f"Worst sampled value: {max(maxima):.2e}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8.5,
            color="#435866", bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75})


def draw_density_pair(fig, reports: dict, left: list[float], right: list[float]) -> None:
    ax_dir = fig.add_axes(left)
    ax_xy = fig.add_axes(right)
    density_panel(ax_dir, reports, "directionUnitVector",
                  "Mapped line / pole: unit-vector error")
    density_panel(ax_xy, reports, "equalAreaUnitHorizon",
                  "Geographic equal-area XY error")


def save_figure(fig, output: Path) -> tuple[Path, Path]:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, facecolor="white")
    pdf = output.with_suffix(".pdf")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    return output, pdf


def make_composite(examples: list, reports: dict, output: Path) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(16.4, 11.1), facecolor="white")
    fig.text(0.045, 0.968, "Geographic rotation: eight routes to the same fabric",
             fontsize=19, fontweight="bold", color="#213748", ha="left", va="top")
    fig.text(0.045, 0.927,
             "All stereonets are lower hemisphere, north up and east right. "
             "Each blue trace and symbol overlays eight plotting conventions and specimen-axis choices.",
             fontsize=10.5, color="#3d5667", ha="left", va="top")

    net_positions = ([0.105, 0.475, 0.31, 0.365], [0.585, 0.475, 0.31, 0.365])
    titles = ("A  Field case from CSV", "B  Synthetic oblique grid case")
    for (field, cases, difference), pos, title in zip(examples, net_positions, titles):
        center = pos[0] + pos[2] / 2
        fig.text(center, 0.874, title, ha="center", va="bottom", fontsize=13,
                 fontweight="bold", color="#213748")
        detail = (f"{field['name']}  ·  strike {field['strike']:g}°, dip {field['dip']:g}° "
                  f"{field['dipDir']}  ·  trend {field['trend']:.1f}°, plunge {field['plunge']:.1f}°")
        fig.text(center, 0.850, detail, ha="center", va="bottom", fontsize=9.2,
                 color="#455d6d")
        ax = fig.add_axes(pos)
        overlay_stereonet(ax, cases)
        fig.text(center, 0.470, f"8 overlays · greatest plotted-coordinate difference {difference:.1e}",
                 ha="center", va="top", fontsize=9, color="#455d6d")

    fig.text(0.5, 0.431,
             "Blue line: foliation great circle       ■  downward foliation pole       "
             "●  lineation (on the plane)",
             ha="center", va="center", fontsize=10.5, color="#344b5d")
    fig.text(0.045, 0.397, "C  Residual density across the full matrix audits",
             fontsize=13, fontweight="bold", color="#213748", va="center")
    fig.text(0.045, 0.375,
             f"MATLAB: {reports['MATLAB MTEX 6.1']['evaluations']:,} mappings  ·  "
             f"Python: {reports['MTEX Python']['total_mapping_cases']:,} mappings  ·  "
             "separate deterministic grids, each spanning eight conventions",
             fontsize=9.7, color="#455d6d", va="center")
    draw_density_pair(fig, reports, [0.065, 0.115, 0.405, 0.22],
                      [0.55, 0.115, 0.405, 0.22])
    fig.text(0.045, 0.065,
             "Each distribution contains every mapping's larger line/pole residual. "
             "XY coordinates use a unit-radius equal-area horizon in both audits. "
             "Values ≤10⁻¹⁷ share the first bin.",
             fontsize=8.8, color="#536a78", va="center")
    fig.text(0.045, 0.039,
             "Default CSV case is synthetic; all specimen axes are illustrative. "
             "Exactly vertical oblique trend-only cases are underdetermined and excluded from the grid.",
             fontsize=8.8, color="#536a78", va="center")
    return save_figure(fig, output)


def make_density_only(reports: dict, output: Path) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(12.5, 5.8), facecolor="white")
    fig.text(0.05, 0.95, "Fabrica rotation residuals across the spherical grid",
             fontsize=17, fontweight="bold", color="#213748", va="top")
    fig.text(0.05, 0.88,
             f"MATLAB {reports['MATLAB MTEX 6.1']['evaluations']:,} mappings  ·  "
             f"Python {reports['MTEX Python']['total_mapping_cases']:,} mappings  ·  "
             "eight conventions in each audit",
             fontsize=10, color="#455d6d", va="top")
    draw_density_pair(fig, reports, [0.07, 0.20, 0.405, 0.56],
                      [0.56, 0.20, 0.405, 0.56])
    fig.text(0.05, 0.075,
             "Fractions are normalized separately; the grids differ. Unit-vector and unit-horizon XY "
             "residuals are dimensionless. Values ≤10⁻¹⁷ share the first bin.",
             fontsize=8.9, color="#536a78", va="center")
    return save_figure(fig, output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--sample", default=DEFAULT_SAMPLE)
    parser.add_argument("--matlab-report", type=Path, default=MATLAB_REPORT)
    parser.add_argument("--python-report", type=Path, default=PYTHON_REPORT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.csv is not None and args.output is None:
        parser.error("--output is required with --csv to preserve the bundled synthetic figures")
    args.csv = args.csv or DEFAULT_CSV
    args.output = args.output or DEFAULT_OUTPUT
    if args.output.suffix.lower() != ".png":
        parser.error("--output must be a .png path")
    reports = load_reports(args.matlab_report, args.python_report)
    examples = plotted_examples(args.csv, args.sample)
    figure_outputs = make_composite(examples, reports, args.output)
    density_outputs = make_density_only(reports, args.output.with_name("geographic_matrix_density.png"))
    for path in (*figure_outputs, *density_outputs):
        print(path)


if __name__ == "__main__":
    main()
