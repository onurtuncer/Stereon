# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Renders the fitting and curve-property figures of the manual
(docs/fitting.rst, docs/properties.rst).

Every curve and number comes from Stereon itself: the program
stereon_figure_data (docs/figures/FigureData.cpp, built with the tests)
prints them as JSON, and this script draws them. Run it after a change to
the fitting or property code, from the repository root:

    python docs/figures/make_fitting_figures.py build/gcc16-debug/tests/stereon_figure_data

or with a JSON file written by that program in place of the executable. The
SVG files are written next to this script and committed, so the manual
builds without a compiler. Requires matplotlib (pip install matplotlib).

Colours follow the manual's chart palette: the first three categorical
slots, validated for colour-vision deficiency and contrast, with neutral ink
for data and construction lines; every series is also named on the figure,
so identity never rests on colour alone.
"""

import json
import os
import subprocess
import sys

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8985"
GRID = "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]  # blue, orange, aqua
BLUE_RAMP = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"]

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK_SECONDARY,
    "axes.edgecolor": INK_MUTED,
    "axes.linewidth": 0.8,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "svg.fonttype": "none",  # text stays text, searchable and small
    "svg.hashsalt": "stereon",  # stable element ids, so reruns diff cleanly
})


def load(source):
    if source.endswith(".json"):
        with open(source, encoding="utf-8") as f:
            return json.load(f)
    output = subprocess.run([source], check=True, capture_output=True, text=True).stdout
    return json.loads(output)


def xy(points):
    return [p[0] for p in points], [p[1] for p in points]


def save(fig, name):
    path = os.path.join(HERE, name)
    fig.savefig(path, format="svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    print("wrote", os.path.relpath(path))


def data_points(ax, points, label="data points", size=22, color=INK):
    x, y = xy(points)
    ax.scatter(x, y, s=size, color=color, zorder=5, label=label, edgecolors=SURFACE,
               linewidths=1.0)


def parametrization(data):
    """Three interpolants of the same points, and where each puts the points."""
    d = data["parametrization"]
    fig, (top, strip) = plt.subplots(
        2, 1, figsize=(7.0, 4.6), gridspec_kw={"height_ratios": [3.2, 1.0]})
    styles = [("uniform", "uniform", "-"), ("chord", "chord length", "--"),
              ("centripetal", "centripetal", "-.")]
    for (key, name, dash), color in zip(styles, SERIES):
        x, y = xy(d[key]["curve"])
        top.plot(x, y, dash, color=color, linewidth=2.0, label=name)
    data_points(top, d["points"])
    top.set_aspect("equal", adjustable="datalim")
    top.set_title("Cubic interpolants of the same eight points", loc="left", fontsize=10)
    top.legend(loc="center", bbox_to_anchor=(0.62, 0.42), ncols=2)
    top.set_xlabel("x (mm)")
    top.set_ylabel("y (mm)")
    for row, ((key, name, _), color) in enumerate(zip(styles, SERIES)):
        t = d[key]["parameters"]
        strip.scatter(t, [row] * len(t), s=26, color=color, zorder=5, edgecolors=SURFACE)
        strip.plot([0, 1], [row, row], color=GRID, linewidth=1.0, zorder=1)
        strip.text(-0.02, row, name, ha="right", va="center", color=INK_SECONDARY)
    strip.set_ylim(-0.6, 2.6)
    strip.set_yticks([])
    strip.set_xlim(-0.01, 1.01)
    strip.set_xlabel("parameter t_k of each point")
    strip.grid(False)
    strip.spines["left"].set_visible(False)
    fig.tight_layout()
    save(fig, "fitting-parametrization.svg")


def knot_averaging(data):
    """Chord-length parameters and the knots averaged from them (eq. 9.8)."""
    d = data["parametrization"]["chord"]
    t = d["parameters"]
    knots = d["knots"]
    degree = 3
    interior = knots[degree + 1:len(knots) - degree - 1]
    fig, ax = plt.subplots(figsize=(7.0, 1.9))
    ax.grid(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([0, 1])
    ax.plot([0, 1], [1, 1], color=INK_MUTED, linewidth=1.0)
    ax.plot([0, 1], [0, 0], color=INK_MUTED, linewidth=1.0)
    for j, u in enumerate(interior, start=1):
        for i in range(j, j + degree):
            ax.plot([t[i], u], [1, 0], color=SERIES[0], linewidth=0.8, alpha=0.45)
    ax.scatter(t, [1] * len(t), s=30, color=INK, zorder=5, edgecolors=SURFACE)
    ax.scatter(interior, [0] * len(interior), s=46, marker="D", color=SERIES[0], zorder=5,
               edgecolors=SURFACE)
    ax.scatter([0, 1], [0, 0], s=46, marker="D", color=SERIES[0], zorder=5, edgecolors=SURFACE)
    ax.text(-0.02, 1, "parameters t_k", ha="right", va="center", color=INK_SECONDARY)
    ax.text(-0.02, 0, "knots u_j", ha="right", va="center", color=INK_SECONDARY)
    ax.text(1.0, 0.5, "each interior knot is the mean of\nthe three parameters joined to it",
            ha="right", va="center", color=INK_SECONDARY, fontsize=8)
    ax.set_ylim(-0.3, 1.3)
    fig.tight_layout()
    save(fig, "fitting-knot-averaging.svg")


def least_squares(data):
    """Least-squares fits with more and more poles, and the interpolant."""
    d = data["leastSquares"]
    fits = d["fits"]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 4.8), sharex=True, sharey=True)
    for ax, fit in zip(axes.flat, fits):
        data_points(ax, d["points"], size=10, color=INK_MUTED)
        px, py = xy(fit["poles"])
        ax.plot(px, py, ":", color=INK_SECONDARY, linewidth=1.0, label="control polygon")
        ax.scatter(px, py, s=18, marker="s", color=INK_SECONDARY, zorder=4)
        x, y = xy(fit["curve"])
        ax.plot(x, y, color=SERIES[0], linewidth=2.0, label="least-squares fit")
        ax.set_title(f"{fit['poleCount']} poles, max deviation {fit['deviation']:.3f} mm",
                     loc="left", fontsize=9)
    ax = axes.flat[3]
    data_points(ax, d["points"], size=10, color=INK_MUTED)
    x, y = xy(d["interpolant"])
    ax.plot(x, y, color=SERIES[1], linewidth=1.4, label="interpolant")
    ax.set_title("interpolant: 60 poles, through every point", loc="left", fontsize=9)
    axes.flat[0].legend(loc="lower left", fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel("x (mm)")
    for ax in axes[:, 0]:
        ax.set_ylabel("y (mm)")
    fig.tight_layout()
    save(fig, "fitting-least-squares.svg")


def approximation(data):
    """Deviation against pole count, the budget, and what Approximate chose."""
    d = data["approximation"]
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    counts = d["counts"]
    deviations = d["deviations"]
    ax.plot(counts, deviations, color=SERIES[0], linewidth=2.0, marker="o", markersize=4,
            markeredgecolor=SURFACE, label="max deviation of the fit")
    ax.axhline(d["budget"], color=INK_SECONDARY, linewidth=1.0, linestyle="--")
    ax.text(counts[-1], d["budget"], "budget", ha="right", va="bottom", color=INK_SECONDARY)
    ax.scatter([d["chosen"]], [d["chosenDeviation"]], s=90, color=SERIES[1], zorder=6,
               edgecolors=SURFACE, linewidths=1.5)
    ax.annotate(f"Approximate chooses {d['chosen']} poles",
                xy=(d["chosen"], d["chosenDeviation"]),
                xytext=(d["chosen"] + 8, d["chosenDeviation"] + 0.12), color=INK,
                arrowprops={"arrowstyle": "-", "color": INK_SECONDARY, "linewidth": 0.8})
    ax.set_xlabel("pole count (cubic least squares)")
    ax.set_ylabel("max deviation at the points (mm)")
    ax.set_ylim(bottom=0)
    ax.set_title("Deviation is not monotone in the pole count", loc="left", fontsize=10)
    fig.tight_layout()
    save(fig, "fitting-approximation.svg")


def collocation(data):
    """The banded collocation matrix N_i(t_k) of a cubic interpolation."""
    rows = data["collocation"]["rows"]
    from matplotlib.colors import ListedColormap, BoundaryNorm

    bounds = [-1e-9, 1e-12, 0.1, 0.25, 0.4, 0.55, 0.7, 1.0001]
    cmap = ListedColormap(BLUE_RAMP)
    norm = BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    image = ax.imshow(rows, cmap=cmap, norm=norm, interpolation="nearest")
    ax.grid(False)
    n = len(rows)
    ax.set_xticks(range(0, n, 2))
    ax.set_yticks(range(0, n, 2))
    ax.set_xticks([x - 0.5 for x in range(n + 1)], minor=True)
    ax.set_yticks([y - 0.5 for y in range(n + 1)], minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=1.0)
    ax.tick_params(which="minor", length=0)
    ax.set_xlabel("pole i")
    ax.set_ylabel("point k")
    ax.set_title("N_i(t_k): at most four non-zeros per row", loc="left", fontsize=10)
    bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, ticks=[0, 0.25, 0.55, 1.0])
    bar.outline.set_visible(False)
    fig.tight_layout()
    save(fig, "fitting-collocation.svg")


def surface(data):
    """A grid on a hemisphere, poles included, and the interpolating surface."""
    d = data["surface"]
    fig = plt.figure(figsize=(5.2, 4.4))
    ax = fig.add_subplot(projection="3d")
    mesh = d["mesh"]
    for row in mesh:
        x, y, z = zip(*row)
        ax.plot(x, y, z, color=SERIES[0], linewidth=0.7)
    for j in range(len(mesh[0])):
        x, y, z = zip(*[row[j] for row in mesh])
        ax.plot(x, y, z, color=SERIES[0], linewidth=0.7)
    x, y, z = zip(*d["grid"])
    ax.scatter(x, y, z, s=14, color=INK, depthshade=False)
    ax.set_title(f"7 x 9 grid interpolated (max deviation {d['deviation']:.0e} mm)",
                 loc="left", fontsize=10)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-60)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRID)
        axis._axinfo["grid"]["color"] = GRID
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_zticks([-1, 0, 1])
    save(fig, "fitting-surface.svg")


def periodic(data):
    """A closed profile fitted periodic and open, and periodic least squares."""
    d = data["periodic"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.4, 3.9))
    x, y = xy(d["open"])
    left.plot(x, y, "--", color=SERIES[1], linewidth=1.6, label="open cubic, first point repeated")
    x, y = xy(d["closed"])
    left.plot(x, y, color=SERIES[0], linewidth=2.0, label="periodic cubic")
    data_points(left, d["points"])
    seam = d["points"][0]
    left.annotate("seam: first and last point", xy=(seam[0], seam[1]),
                  xytext=(seam[0] + 0.55, seam[1] + 0.55), color=INK, fontsize=8,
                  arrowprops={"arrowstyle": "-", "color": INK_SECONDARY, "linewidth": 0.8})
    left.set_aspect("equal", adjustable="datalim")
    left.set_title("Twelve points on a cam profile", loc="left", fontsize=10)
    left.set_xlabel("x (mm)")
    left.set_ylabel("y (mm)")
    left.legend(loc="lower left", fontsize=8)
    data_points(right, d["noisy"], size=10, color=INK_MUTED)
    px, py = xy(d["poles"])
    right.plot(px, py, ":", color=INK_SECONDARY, linewidth=1.0, label="control polygon")
    right.scatter(px, py, s=18, marker="s", color=INK_SECONDARY, zorder=4)
    x, y = xy(d["interpolant"])
    right.plot(x, y, color=SERIES[1], linewidth=1.0, label="periodic interpolant")
    x, y = xy(d["leastSquares"])
    right.plot(x, y, color=SERIES[0], linewidth=2.0,
               label=f"periodic least squares, {d['poleCount']} poles")
    right.set_aspect("equal", adjustable="datalim")
    right.set_title(f"Fifty noisy samples, max deviation {d['deviation']:.3f} mm",
                    loc="left", fontsize=10)
    right.set_xlabel("x (mm)")
    right.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    save(fig, "fitting-periodic.svg")


def curvature_comb(data):
    """The curvature comb: teeth along the principal normal, scaled by curvature."""
    d = data["properties"]
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    for base, tip in zip(d["curve"], d["comb"]):
        ax.plot([base[0], tip[0]], [base[1], tip[1]], color=SERIES[0], linewidth=0.6,
                alpha=0.7)
    x, y = xy(d["comb"])
    ax.plot(x, y, color=SERIES[0], linewidth=1.2, label="curvature comb")
    x, y = xy(d["curve"])
    ax.plot(x, y, color=INK, linewidth=2.0, label="cubic interpolant")
    data_points(ax, d["points"], color=INK_SECONDARY)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_title("Curvature comb of a cubic through eight points (longest tooth 1.5 mm)",
                 loc="left", fontsize=10)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.legend(loc="lower right")
    fig.tight_layout()
    save(fig, "properties-comb.svg")


def arc_length(data):
    """Stations at equal parameter steps and at equal arc length, and s(t)."""
    d = data["properties"]
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(7.4, 3.4), gridspec_kw={"width_ratios": [1.7, 1.0]})
    x, y = xy(d["curve"])
    left.plot(x, y, color=INK_MUTED, linewidth=1.4, zorder=1)
    x, y = xy(d["by_parameter"])
    left.scatter(x, y, s=34, color=SERIES[1], marker="s", zorder=4, edgecolors=SURFACE,
                 label="equal steps of t")
    x, y = xy(d["by_length"])
    left.scatter(x, y, s=30, color=SERIES[0], zorder=5, edgecolors=SURFACE,
                 label="equal steps of arc length")
    left.set_aspect("equal", adjustable="datalim")
    left.set_title("Thirteen stations, two ways", loc="left", fontsize=10)
    left.set_xlabel("x (mm)")
    left.set_ylabel("y (mm)")
    left.legend(loc="lower right")
    lengths = d["lengths_at_parameter"]
    t = [i / (len(lengths) - 1) for i in range(len(lengths))]
    right.plot(t, lengths, "-", color=SERIES[1], linewidth=2.0, marker="s", markersize=4,
               label="s(t)")
    right.plot([0, 1], [0, d["length"]], "--", color=INK_MUTED, linewidth=1.0,
               label="constant speed")
    right.set_title(f"Arc length, total {d['length']:.4f} mm", loc="left", fontsize=10)
    right.set_xlabel("parameter t")
    right.set_ylabel("s (mm)")
    right.legend(loc="upper left")
    fig.tight_layout()
    save(fig, "properties-arclength.svg")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    data = load(sys.argv[1])
    parametrization(data)
    knot_averaging(data)
    least_squares(data)
    approximation(data)
    collocation(data)
    surface(data)
    periodic(data)
    curvature_comb(data)
    arc_length(data)


if __name__ == "__main__":
    main()
