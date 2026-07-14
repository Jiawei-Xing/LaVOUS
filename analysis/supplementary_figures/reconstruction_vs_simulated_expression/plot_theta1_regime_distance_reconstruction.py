#!/usr/bin/env python
"""Reconstruction diagnostics by regime and root distance across theta settings."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
GENE_ID = 500
TREE_PATH = REPO_ROOT / "lineage_simulation" / "tree.nwk"

SCENARIOS = [
    {"scenario": "Base", "theta1": 1, "title": r"$\theta_1=1$"},
    {"scenario": "Theta1", "theta1": 3, "title": r"$\theta_1=3$"},
    {"scenario": "t5r5", "theta1": 5, "title": r"$\theta_1=5$"},
    {"scenario": "t7r5", "theta1": 7, "title": r"$\theta_1=7$"},
]

REGIME_COLORS = {
    0: "#9a9a9a",
    1: "#c43c39",
}
DISTANCE_BINS = [-0.1, 4, 8, 12, 16, 20, 24, 28, 32.1]
DISTANCE_BIN_LABELS = [
    "0-4",
    "4-8",
    "8-12",
    "12-16",
    "16-20",
    "20-24",
    "24-28",
    "28-32",
]


def tree_depths(tree_path: Path) -> pd.DataFrame:
    tree = Phylo.read(tree_path, "newick")
    rows = [
        {"node_name": clade.name, "root_distance": float(distance)}
        for clade, distance in tree.depths().items()
        if clade.name
    ]
    return pd.DataFrame(rows)


def load_scenario(config: dict[str, object], depths: pd.DataFrame) -> pd.DataFrame:
    scenario = str(config["scenario"])
    theta1 = int(config["theta1"])
    sim_path = (
        REPO_ROOT
        / "expression_simulation"
        / "simulation"
        / f"sim_history_{scenario}.tsv"
    )
    reconst_path = (
        REPO_ROOT
        / "expression_simulation"
        / "reconst"
        / f"history_{scenario}_gene{GENE_ID}.tsv"
    )

    sim = pd.read_csv(sim_path, sep="	")
    reconst = pd.read_csv(reconst_path, sep="	")
    merged = sim.merge(
        reconst[["node_name", "regime", "infer_mu", "infer_var"]],
        on="node_name",
        how="inner",
        validate="one_to_one",
    ).merge(depths, on="node_name", how="inner", validate="one_to_one")

    if len(merged) != len(sim) or len(merged) != len(reconst):
        raise ValueError(
            f"{scenario}: matched {len(merged)} nodes, "
            f"but simulation has {len(sim)} and reconstruction has {len(reconst)}"
        )

    sim_root = float(sim["sim_expr"].iloc[0])
    merged.insert(0, "scenario", scenario)
    merged.insert(1, "theta1", theta1)
    merged.insert(2, "theta_label", str(config["title"]))
    merged["gene_id"] = GENE_ID
    merged["sim_root"] = sim_root
    merged["sim_centered"] = merged["sim_expr"].astype(float) - sim_root
    merged["reconstructed_centered"] = merged["infer_mu"].astype(float) - sim_root
    merged["residual"] = merged["infer_mu"].astype(float) - merged["sim_expr"].astype(
        float
    )
    merged["abs_error"] = merged["residual"].abs()
    merged["regime"] = merged["regime"].astype(int)
    merged["regime_label"] = np.where(merged["regime"] == 0, "Regime 0", "Regime 1")
    merged["tree_file"] = str(TREE_PATH)
    merged["sim_history_file"] = str(sim_path)
    merged["reconstruction_file"] = str(reconst_path)
    return merged


def pearsonr(x: pd.Series, y: pd.Series) -> float:
    if len(x) < 2:
        return np.nan
    return float(
        np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1]
    )


def summarize(source_data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scenario, theta1, regime), group in source_data.groupby(
        ["scenario", "theta1", "regime"], sort=False
    ):
        for node_class, sub in [
            ("all nodes", group),
            ("internal nodes", group[~group["is_leaf"].astype(bool)]),
            ("leaf nodes", group[group["is_leaf"].astype(bool)]),
        ]:
            residual = sub["residual"].astype(float)
            rows.append(
                {
                    "scenario": scenario,
                    "theta1": theta1,
                    "gene_id": GENE_ID,
                    "regime": regime,
                    "node_class": node_class,
                    "n_nodes": len(sub),
                    "pearson_r": pearsonr(
                        sub["sim_centered"], sub["reconstructed_centered"]
                    ),
                    "mae": float(np.mean(np.abs(residual))),
                    "rmse": float(np.sqrt(np.mean(residual**2))),
                    "bias": float(np.mean(residual)),
                    "mean_sim_expr": float(np.mean(sub["sim_expr"])),
                    "mean_reconstructed_expr": float(np.mean(sub["infer_mu"])),
                    "median_infer_var": float(np.median(sub["infer_var"])),
                    "mean_infer_var": float(np.mean(sub["infer_var"])),
                    "min_root_distance": float(np.min(sub["root_distance"])),
                    "max_root_distance": float(np.max(sub["root_distance"])),
                }
            )
    return pd.DataFrame(rows)


def binned_error_summary(source_data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scenario, theta1), group in source_data.groupby(
        ["scenario", "theta1"], sort=False
    ):
        binned = group.copy()
        binned["distance_bin"] = pd.cut(
            binned["root_distance"],
            bins=DISTANCE_BINS,
            labels=DISTANCE_BIN_LABELS,
            include_lowest=True,
        )
        for distance_bin, sub in binned.groupby("distance_bin", observed=True):
            rows.append(
                {
                    "scenario": scenario,
                    "theta1": theta1,
                    "gene_id": GENE_ID,
                    "distance_bin": str(distance_bin),
                    "n_nodes": len(sub),
                    "median_root_distance": float(np.median(sub["root_distance"])),
                    "median_abs_error": float(np.median(sub["abs_error"])),
                    "q25_abs_error": float(np.quantile(sub["abs_error"], 0.25)),
                    "q75_abs_error": float(np.quantile(sub["abs_error"], 0.75)),
                }
            )
    return pd.DataFrame(rows)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.18,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
    )


def set_equal_limits(ax: plt.Axes, x: pd.Series, y: pd.Series) -> None:
    values = np.concatenate([np.asarray(x, dtype=float), np.asarray(y, dtype=float)])
    lo = float(np.nanpercentile(values, 0.5))
    hi = float(np.nanpercentile(values, 99.5))
    pad = 0.08 * max(hi - lo, 1.0)
    lo -= pad
    hi += pad
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ticks = ax.get_yticks()
    ticks = ticks[(ticks >= lo) & (ticks <= hi)]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.plot([lo, hi], [lo, hi], color="#303030", lw=0.8, ls="--", zorder=1)


def set_pm6_ticks(ax: plt.Axes) -> None:
    ticks = np.array([-6, -4, -2, 0, 2, 4, 6], dtype=float)
    lo = min(ax.get_xlim()[0], ax.get_ylim()[0], -6.2)
    hi = max(ax.get_xlim()[1], ax.get_ylim()[1], 6.2)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.plot([lo, hi], [lo, hi], color="#303030", lw=0.8, ls="--", zorder=1)


def set_y_max_tick(ax: plt.Axes, max_tick: float) -> None:
    lo = min(ax.get_ylim()[0], 0.0)
    hi = max(ax.get_ylim()[1], max_tick + 0.05 * max(max_tick - lo, 1.0))
    ax.set_ylim(lo, hi)
    ax.set_yticks(np.arange(0, max_tick + 0.001, 1.0))


def scatter_by_regime(
    ax: plt.Axes,
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    *,
    s: float,
    alpha: float,
    legend: bool = False,
) -> None:
    for regime in sorted(data["regime"].unique()):
        sub = data[data["regime"] == regime]
        ax.scatter(
            sub[x_col],
            sub[y_col],
            s=s,
            color=REGIME_COLORS.get(int(regime), "#4d4d4d"),
            alpha=alpha,
            edgecolors="none",
            rasterized=True,
            label=f"Regime {int(regime)}" if legend else "_nolegend_",
            zorder=3 if int(regime) != 0 else 2,
        )
    if legend:
        ax.legend(frameon=False, loc="lower right", handletextpad=0.2, borderpad=0.1)


def add_binned_median_line(
    ax: plt.Axes, data: pd.DataFrame, *, legend: bool = False
) -> None:
    binned = binned_error_summary(data)
    if binned.empty:
        return
    (line,) = ax.plot(
        binned["median_root_distance"],
        binned["median_abs_error"],
        color="#202020",
        lw=1.2,
        marker="o",
        markersize=3.2,
        markerfacecolor="white",
        markeredgecolor="#202020",
        markeredgewidth=0.8,
        label="Binned median",
        zorder=5,
    )
    if legend:
        ax.legend(
            handles=[line],
            labels=["Binned median"],
            frameon=False,
            loc="upper left",
            handletextpad=0.3,
            borderpad=0.1,
        )


def add_accuracy_label(ax: plt.Axes, data: pd.DataFrame) -> None:
    corr = pearsonr(data["sim_centered"], data["reconstructed_centered"])
    rmse = float(np.sqrt(np.mean(np.square(data["residual"].astype(float)))))
    ax.text(
        0.05,
        0.95,
        f"$r$={corr:.2f}\nRMSE={rmse:.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        linespacing=1.15,
    )


def plot_figure(source_data: pd.DataFrame, out_prefix: Path) -> None:
    mpl.rcParams.update(
        {
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.7,
        }
    )

    fig, axes = plt.subplots(
        len(SCENARIOS), 3, figsize=(7.6, 9.2), constrained_layout=False
    )
    letters = list("ABCDEFGHIJKL")

    for row_idx, config in enumerate(SCENARIOS):
        scenario = str(config["scenario"])
        data = source_data[source_data["scenario"] == scenario].copy()
        row_axes = axes[row_idx, :]
        theta_title = str(config["title"])

        ax = row_axes[0]
        scatter_by_regime(
            ax,
            data,
            "sim_centered",
            "reconstructed_centered",
            s=8,
            alpha=0.62,
            legend=row_idx == 0,
        )
        set_equal_limits(ax, data["sim_centered"], data["reconstructed_centered"])
        set_pm6_ticks(ax)
        add_accuracy_label(ax, data)
        ax.set_title(theta_title)
        ax.set_ylabel("Reconstructed expression")
        ax.set_xlabel("Simulated expression")

        ax = row_axes[1]
        scatter_by_regime(
            ax,
            data,
            "root_distance",
            "abs_error",
            s=7,
            alpha=0.50,
        )
        add_binned_median_line(ax, data, legend=row_idx == 0)
        ax.set_title("Reconstruction error")
        ax.set_ylabel("Absolute error")
        ax.set_xlabel("Distance from root")
        if row_idx == 3:
            set_y_max_tick(ax, 5.0)

        ax = row_axes[2]
        scatter_by_regime(
            ax,
            data,
            "root_distance",
            "infer_var",
            s=7,
            alpha=0.50,
        )
        ax.set_title("Reconstruction variance")
        ax.set_ylabel("Posterior variance")
        ax.set_xlabel("Distance from root")
        if row_idx == 0:
            set_y_max_tick(ax, 6.0)
        elif row_idx in (2, 3):
            set_y_max_tick(ax, 7.0)

    for ax, letter in zip(axes.ravel(), letters):
        add_panel_label(ax, letter)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(length=2.5, width=0.7)

    fig.subplots_adjust(
        left=0.075, right=0.985, top=0.935, bottom=0.065, wspace=0.36, hspace=0.60
    )
    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".png"), dpi=450)
    plt.close(fig)


def write_caption(
    outdir: Path, out_prefix: Path, source_data: pd.DataFrame, summary: pd.DataFrame
) -> None:
    stats = []
    for config in SCENARIOS:
        scenario = str(config["scenario"])
        theta1 = int(config["theta1"])
        sub = source_data[source_data["scenario"] == scenario]
        stats.append(
            f"theta1={theta1}: r="
            f"{pearsonr(sub['sim_centered'], sub['reconstructed_centered']):.2f}, "
            f"MAE={sub['abs_error'].mean():.2f}"
        )
    source_files = source_data[
        ["tree_file", "sim_history_file", "reconstruction_file"]
    ].drop_duplicates()
    file_lines = "\n".join(
        [
            f"- `{row.tree_file}`, `{row.sim_history_file}`, and `{row.reconstruction_file}`"
            for row in source_files.itertuples(index=False)
        ]
    )
    caption = f"""# Supplementary figure: reconstruction diagnostics by theta on the original r=5 tree simulations

Panels A-C show the theta1=1, r=5 Gene 500 simulation on the original tree; panels D-F show theta1=3, r=5; panels G-I show theta1=5, r=5; panels J-L show theta1=7, r=5. In each row, the first panel compares LaVOUS posterior mean reconstructed expression with simulated latent expression at matched tree nodes, the second panel shows absolute reconstruction error versus distance from the root, and the third panel shows posterior reconstruction variance versus distance from the root. Black open circles and connecting lines in the error panels show median absolute error within fixed 4-unit root-distance bins. Expression values in the first column are centered by the simulated root expression for each scenario. Points are colored by reconstructed branch regime: gray for regime 0 and red for regime 1.

Summary across all nodes: {'; '.join(stats)}.

Source files:
{file_lines}

Generated files:
- `{out_prefix.with_suffix('.pdf').name}`
- `{out_prefix.with_suffix('.png').name}`
- `theta1_regime_distance_source_data.tsv`
- `theta1_regime_distance_summary.tsv`
- `theta1_regime_distance_error_bins.tsv`
"""
    (outdir / "caption_theta1_regime_distance.md").write_text(caption)


def main() -> None:
    THIS_DIR.mkdir(parents=True, exist_ok=True)
    depths = tree_depths(TREE_PATH)
    source_data = pd.concat(
        [load_scenario(config, depths) for config in SCENARIOS], ignore_index=True
    )
    summary = summarize(source_data)

    source_data.to_csv(
        THIS_DIR / "theta1_regime_distance_source_data.tsv", sep="	", index=False
    )
    summary.to_csv(
        THIS_DIR / "theta1_regime_distance_summary.tsv", sep="	", index=False
    )
    binned_error_summary(source_data).to_csv(
        THIS_DIR / "theta1_regime_distance_error_bins.tsv", sep="	", index=False
    )

    out_prefix = THIS_DIR / "lavous_theta1_regime_distance_reconstruction"
    plot_figure(source_data, out_prefix)
    write_caption(THIS_DIR, out_prefix, source_data, summary)


if __name__ == "__main__":
    main()
