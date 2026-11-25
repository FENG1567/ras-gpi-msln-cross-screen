"""Render the four main figures from the released figure-level summaries."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


# Muted, high-contrast colors remain distinguishable when printed or viewed in
# grayscale. Marker shape also identifies each screen in Figure 1.
BLUE = "#176B87"
ORANGE = "#B5663C"
SLATE = "#455A64"
GRAY = "#7B858C"
LIGHT_GRAY = "#AAB2B8"
DARK = "#1D303D"
GRID = "#E4E9EC"

WIDTH_MM = 170


def configure_plotting(arial_bold):
    if arial_bold is not None:
        font_manager.fontManager.addfont(arial_bold)
        font = font_manager.FontProperties(fname=arial_bold)
        if font.get_name() != "Arial":
            raise ValueError("--arial-bold must point to an installed Arial Bold font")
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 9,
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.labelsize": 9,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "xtick.labelsize": 8.8,
            "ytick.labelsize": 8.8,
            "xtick.major.pad": 2.2,
            "ytick.major.pad": 2.0,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "legend.frameon": False,
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def figure(height_mm):
    return plt.figure(figsize=(WIDTH_MM / 25.4, height_mm / 25.4))


def axes_at(fig, left_mm, bottom_mm, width_mm, height_mm):
    page_h = fig.get_figheight() * 25.4
    return fig.add_axes([left_mm / WIDTH_MM, bottom_mm / page_h,
                         width_mm / WIDTH_MM, height_mm / page_h])


def panel(fig, x_mm, y_mm, letter, title):
    y = y_mm / (fig.get_figheight() * 25.4)
    fig.text(x_mm / WIDTH_MM, y, letter, ha="left", va="center",
             fontsize=18, fontweight="bold", color=DARK)
    fig.text((x_mm + 8) / WIDTH_MM, y, title, ha="left", va="center",
             fontsize=10.5, fontweight="bold", color=DARK)


def save_plot(fig, base, previews):
    base.parent.mkdir(parents=True, exist_ok=True)
    previews.mkdir(parents=True, exist_ok=True)
    tiff = base.parent / "TIFF"
    tiff.mkdir(parents=True, exist_ok=True)
    editable = base.parent / "SVG"
    editable.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(editable / f"{base.name}.svg")
    fig.savefig(previews / f"{base.name}.png", dpi=300)
    fig.savefig(tiff / f"{base.name}.tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def plot_fig1(immune, candidates, out, preview):
    fig = figure(76)
    a = axes_at(fig, 30, 21.5, 34, 41.5)
    b = axes_at(fig, 104, 21.5, 62, 41.5)

    # Separate the GPI-panel fraction from the much smaller all-gene baseline
    # so the baseline count and denominator remain legible on the same scale.
    row_y = []
    row_labels = []
    row_values = []
    row_colors = []
    row_texts = []
    for screen_index, row in enumerate(immune[:2]):
        for offset, numerator, denominator, label, color in (
            (0.0, row["hits"], row["gpi"], "GPI panel", BLUE),
            (
                0.82,
                row["background_hits"],
                row["n"],
                "All scored genes",
                GRAY,
            ),
        ):
            y = screen_index * 2.2 + offset
            pct = 100 * numerator / denominator
            pct_digits = 1 if label == "GPI panel" else 2
            row_y.append(y)
            row_labels.append(f"{row['label']}\n{label}")
            row_values.append(pct)
            row_colors.append(color)
            row_texts.append(
                f"{numerator:,}/{denominator:,}\n({pct:.{pct_digits}f}%)"
            )

    for y, pct, color, text in zip(row_y, row_values, row_colors, row_texts):
        a.barh(y, pct, height=0.46, color=color, zorder=2)
        a.text(1.055, y, text, transform=a.get_yaxis_transform(),
               va="center", ha="left", fontsize=8.5, linespacing=1.2)
    a.set_yticks(row_y, row_labels)
    a.invert_yaxis()
    a.set_xlim(-1.5, 80)
    a.set_xticks([0, 20, 40, 60, 80])
    a.margins(y=0.13)
    a.set_xlabel("Genes meeting escape criteria (%)")
    a.grid(axis="x", color=GRID, linewidth=0.55, zorder=0)
    a.set_axisbelow(True)
    a.axhline(1.52, color=LIGHT_GRAY, lw=0.65, ls=(0, (2, 2)))
    panel(fig, 3, 70, "A", "MSLN screen enrichment")

    series = [
        ("drug", "AsPC-1", ORANGE, "o", -0.22),
        ("hup", "HuP-T3", BLUE, "^", 0.0),
        ("kp4", "KP4_MSLN", SLATE, "s", 0.22),
    ]
    for key, label, color, marker, shift in series:
        for i, item in enumerate(candidates):
            value = item.get(key)
            if value is None:
                continue
            hit = item.get(key + "_hit", True)
            b.scatter(
                value,
                i + shift,
                s=25,
                marker=marker,
                edgecolor=color,
                facecolor=color if hit else "white",
                linewidth=0.95,
                zorder=3,
            )

    b.set_yticks(range(len(candidates)), [row["gene"] for row in candidates])
    b.invert_yaxis()
    b.set_xlim(0, 4.85)
    b.set_xlabel("Signed screen effect (log2 units)")
    b.axvline(0.5, color=GRAY, ls=(0, (2, 2)), lw=0.8, zorder=1)
    b.grid(axis="x", color=GRID, linewidth=0.55, zorder=0)
    b.set_axisbelow(True)
    b.margins(y=0.14)
    panel(fig, 87, 70, "B", "Cross-screen candidate effects")

    legend_handles = [
        Line2D([], [], marker=marker, linestyle="none", color=color, ms=4.7,
               label=label)
        for _, label, color, marker, _ in series
    ]
    legend_handles.extend(
        [
            Line2D(
                [], [], marker="o", linestyle="none", markeredgecolor=DARK,
                markerfacecolor=DARK, color=DARK, ms=4.5,
                label="Hit",
            ),
            Line2D(
                [], [], marker="o", linestyle="none", markeredgecolor=DARK,
                markerfacecolor="white", color=DARK, ms=4.5,
                label="Non-hit",
            ),
            Line2D(
                [], [], color=GRAY, linestyle=(0, (2, 2)), lw=0.8,
                label="Effect cut-off (0.5)",
            ),
        ]
    )
    fig.legend(
        handles=[legend_handles[i] for i in (0, 3, 1, 4, 2, 5)],
        loc="center",
        bbox_to_anchor=(0.745, 0.079),
        ncol=3,
        prop={"weight": "bold", "size": 8.2},
        handletextpad=0.28,
        columnspacing=0.7,
        labelspacing=0.25,
        borderaxespad=0.0,
    )
    save_plot(fig, out / "Fig1_MSLN_escape", preview)


def intervals(axis, rows, xlabel, annotate_n=True):
    for i, row in enumerate(rows):
        color = BLUE if "GPI" in row["name"] else GRAY
        axis.plot(
            [row["q1"], row["q3"]], [i, i], color=color, lw=2.2,
            solid_capstyle="round", zorder=2,
        )
        axis.plot(
            row["median"], i, marker="o", color=DARK, markeredgecolor="white",
            markeredgewidth=0.45, ms=4.7, zorder=3,
        )
    labels = [row["name"].replace("RPL/RPS prefix", "RPL/RPS genes") for row in rows]
    if annotate_n:
        labels = [f"{name}\n(n={row['n']:,})" for name, row in zip(labels, rows)]
    axis.set_yticks(range(len(rows)), labels)
    axis.invert_yaxis()
    axis.axvline(0, color=LIGHT_GRAY, ls=(0, (2, 2)), lw=0.8, zorder=1)
    axis.grid(axis="x", color=GRID, linewidth=0.55, zorder=0)
    axis.set_axisbelow(True)
    axis.set_xlabel(xlabel)
    axis.margins(x=0.09, y=0.28)


def plot_fig4(drug, depmap, out, preview):
    fig = figure(55)
    a = axes_at(fig, 31, 15, 48, 29)
    b = axes_at(fig, 116, 15, 50, 29)
    intervals(a, drug, "Relative depletion under\nMRTX1133 (log2 ratio)")
    intervals(b, depmap, "Chronos gene effect")
    a.set_xticks([-0.4, 0, 0.4, 0.8])
    b.set_xticks([-2, -1.5, -1, -0.5, 0])
    panel(fig, 3, 49, "A", "AsPC-1 drug screen")
    panel(fig, 86, 49, "B", "Baseline dependency")
    save_plot(fig, out / "Fig4_drug_and_dependency", preview)


def plot_fig3(rows, out, preview):
    fig = figure(65)
    ax = axes_at(fig, 62, 10, 104, 51)
    positions = [0, 1, 2.55, 3.55, 4.55, 5.55, 6.55]
    label_flags = {"Fröse p1112": "†", "U87": "†", "Nalm6 day 1": "*"}
    for index, (row, y) in enumerate(zip(rows, positions)):
        pct = 100 * row["hits"] / row["gpi"]
        color = BLUE if index < 2 else GRAY
        ax.hlines(y, 0, pct, lw=2.2, color=color, zorder=2)
        ax.plot(pct, y, "o", ms=4.3, color=color, zorder=3)
        ax.text(max(pct + 1.45, 1.25), y, f"{row['hits']}/{row['gpi']}",
                va="center", ha="left", fontsize=9)
    ax.set(
        xlabel="Covered GPI genes meeting escape criteria (%)",
        xlim=(-1.5, 80),
        xticks=range(0, 81, 10),
        yticks=positions,
        yticklabels=[
            f"{row['label']}{label_flags.get(row['label'], '')} · {row['antigen']}"
            for row in rows
        ],
    )
    ax.invert_yaxis()
    ax.grid(axis="x", color=GRID, linewidth=0.55, zorder=0)
    ax.set_axisbelow(True)
    ax.axhline(1.78, color=LIGHT_GRAY, lw=0.65, ls=(0, (2, 2)))
    save_plot(fig, out / "Fig3_immune_comparisons", preview)


def plot_fig2(grid, nulls, out, preview):
    fig = figure(76)
    a = axes_at(fig, 24, 18, 47, 45)
    b = axes_at(fig, 117, 18, 49, 45)
    pairs = sorted({(row["emin"], row["pmin"]) for row in grid})
    xs = sorted({row["dmin"] for row in grid})
    lookup = {
        (row["dmin"], row["emin"], row["pmin"]): row["hits"] for row in grid
    }
    counts = [[lookup[(x, *group)] for x in xs] for group in pairs]
    xedges = [i - 0.5 for i in range(len(xs) + 1)]
    yedges = [i - 0.5 for i in range(len(pairs) + 1)]
    vmax = max(map(max, counts)) + 4
    a.pcolormesh(
        xedges,
        yedges,
        counts,
        cmap="Blues",
        vmin=0,
        vmax=vmax,
        shading="flat",
        rasterized=False,
    )
    a.invert_yaxis()
    for i, group in enumerate(counts):
        for j, value in enumerate(group):
            rgb = matplotlib.colormaps["Blues"](value / vmax)[:3]
            linear = [c / 12.92 if c < 0.04045 else ((c + .055) / 1.055) ** 2.4 for c in rgb]
            luminance = sum(c * weight for c, weight in zip(linear, (.2126, .7152, .0722)))
            color = "white" if luminance < .179 else "#10161A"
            a.text(j, i, str(value), va="center", ha="center", color=color, fontsize=10)
    a.set(
        xticks=range(len(xs)),
        xticklabels=[f"{x:g}" for x in xs],
        yticks=range(len(pairs)),
        yticklabels=[f"{e:g} / {p:g}" for e, p in pairs],
        xlabel="Minimum depletion score",
        ylabel="Minimum CAR log2 fold change\n/ source score",
    )
    baseline = (0.5, 0.5, 2.0)
    baseline_cell = (xs.index(baseline[0]), pairs.index((baseline[1], baseline[2])))
    a.add_patch(
        Rectangle(
            (baseline_cell[0] - 0.48, baseline_cell[1] - 0.48),
            0.96,
            0.96,
            fill=False,
            edgecolor=ORANGE,
            linewidth=1.8,
            zorder=4,
        )
    )
    for tick, value in zip(a.get_xticklabels(), xs):
        if value == baseline[0]:
            tick.set_fontweight("bold")
            tick.set_color(ORANGE)
    for tick, pair in zip(a.get_yticklabels(), pairs):
        if pair == baseline[1:]:
            tick.set_fontweight("bold")
            tick.set_color(ORANGE)
    a.grid(False)
    panel(fig, 3, 70, "A", "GPI dual hits")

    for i, row in enumerate(nulls):
        # q1/q3 in this summary map to the source q025/q975 endpoints, so these
        # bars encode central 95% random-set intervals rather than an IQR.
        b.plot(
            [row["q1"], row["q3"]],
            [i, i],
            color=LIGHT_GRAY,
            lw=2.2,
            solid_capstyle="round",
            zorder=1,
        )
        b.plot(row["median"], i, "o", color=GRAY, ms=4.4, zorder=2)
        b.plot(row["observed"], i, "D", color=BLUE, ms=4.4, zorder=3)
    match_names = {
        "CPM + CAR, guide median": "CPM + both CAR ranks\n(guide median)",
        "CAR neighbours, guide median": "CPM + CAR neighbours\n(guide median)",
    }
    row_labels = [
        match_names.get(row["name"], row["name"].replace(", guide median", "\n(guide median)"))
        for row in nulls
    ]
    b.set(
        yticks=range(len(nulls)),
        yticklabels=row_labels,
        xlabel="Mean MRTX1133 depletion rank",
        xlim=(0, 1),
        xticks=[0, 0.25, 0.5, 0.75, 1.0],
        xticklabels=["0", "0.25", "0.50", "0.75", "1.00"],
    )
    b.invert_yaxis()
    b.grid(axis="x", color=GRID, linewidth=0.55, zorder=0)
    b.set_axisbelow(True)
    b.margins(y=0.25)
    for tick in b.get_yticklabels():
        if tick.get_text() == "CPM + both CAR ranks":
            tick.set_color(BLUE)
    panel(fig, 77, 70, "B", "Matched random sets")
    legend_handles = [
        Line2D([], [], marker="D", linestyle="none", color=BLUE, ms=4.4,
               label="GPI mean"),
        Line2D([], [], marker="o", linestyle="none", color=GRAY, ms=4.4,
               label="Matched median"),
        Line2D([], [], color=LIGHT_GRAY, lw=2.2, solid_capstyle="round",
               label="Central 95% interval"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="center right",
        bbox_to_anchor=(0.981, 0.066),
        ncol=3,
        prop={"weight": "bold", "size": 8.2},
        handletextpad=0.45,
        columnspacing=0.8,
        borderaxespad=0.0,
    )
    save_plot(fig, out / "Fig2_overlap_and_rank", preview)


def load_table(path, floats=(), integers=(), booleans=()):
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    for row in rows:
        for key in floats:
            row[key] = float(row[key]) if row[key] else None
        for key in integers:
            row[key] = int(row[key])
        for key in booleans:
            row[key] = row[key].casefold() in {"true", "1"}
    return rows


def default_data_dir():
    script_dir = Path(__file__).resolve().parent
    candidates = (
        script_dir.parent / "data" / "derived",
        script_dir.parent / "readiness_repo" / "data" / "derived",
    )
    return next((path for path in candidates if path.is_dir()), candidates[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=default_data_dir())
    parser.add_argument("--output-dir", type=Path, default=Path("generated/figures"))
    parser.add_argument("--arial-bold", type=Path,
                        help="Use a locally licensed Arial Bold font file to match the manuscript typography")
    args = parser.parse_args()
    data, out = args.data, args.output_dir
    preview = out / "png"
    preview.mkdir(parents=True, exist_ok=True)
    configure_plotting(args.arial_bold)

    immune = load_table(
        data / "immune_screen_summary.tsv",
        integers=("hits", "gpi", "background_hits", "n"),
    )
    candidates = load_table(
        data / "primary_candidate_effects.tsv",
        floats=("drug", "hup", "kp4"),
        booleans=("drug_hit", "hup_hit", "kp4_hit"),
    )
    dist_opts = dict(floats=("q1", "median", "q3"), integers=("n",))
    drug = load_table(data / "drug_distribution.tsv", **dist_opts)
    depmap = load_table(data / "depmap_baseline_summary.tsv", **dist_opts)
    grid = load_table(
        data / "primary_threshold_grid.tsv",
        floats=("dmin", "emin", "pmin"),
        integers=("hits",),
    )
    nulls = load_table(
        data / "matched_rank_summary.tsv",
        floats=("observed", "median", "q1", "q3"),
    )

    plot_fig1(immune, candidates, out, preview)
    plot_fig2(grid, nulls, out, preview)
    plot_fig3(immune, out, preview)
    plot_fig4(drug, depmap, out, preview)
    print("Rendered four 170-mm figures as vector PDF/SVG and 300-dpi PNG/TIFF.")


if __name__ == "__main__":
    main()
