#!/usr/bin/env python
"""Compare LaVOUS reconstructed and simulated expression histories from Fig. 2D."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCENARIOS = [
    ("r0", 1),
    ("t3r0", 3),
    ("t5r0", 5),
    ("t7r0", 7),
]


def load_scenario(simdir: Path, reconstdir: Path, scenario: str, theta1: int) -> pd.DataFrame:
    sim_path = simdir / f"sim_history_{scenario}.tsv"
    reconst_path = reconstdir / f"history_{scenario}_gene500.tsv"

    sim = pd.read_csv(sim_path, sep="\t")
    reconst = pd.read_csv(reconst_path, sep="\t")
    merged = sim.merge(
        reconst[["node_name", "infer_mu", "infer_var"]],
        on="node_name",
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != len(sim) or len(merged) != len(reconst):
        raise ValueError(
            f"{scenario}: matched {len(merged)} nodes, "
            f"but simulation has {len(sim)} and reconstruction has {len(reconst)}"
        )

    sim_root = float(sim["sim_expr"].iloc[0])
    merged.insert(0, "scenario", scenario)
    merged.insert(1, "theta1", theta1)
    merged["sim_root"] = sim_root
    merged["sim_centered"] = merged["sim_expr"].astype(float) - sim_root
    merged["reconstructed_centered"] = merged["infer_mu"].astype(float) - sim_root
    merged["residual"] = merged["reconstructed_centered"] - merged["sim_centered"]
    merged["abs_error"] = merged["residual"].abs()
    merged["sim_history_file"] = str(sim_path)
    merged["reconstruction_file"] = str(reconst_path)
    return merged


def pearsonr(x: pd.Series, y: pd.Series) -> float:
    if len(x) < 2:
        return np.nan
    return float(np.corrcoef(np.asarray(x, dtype=float), np.asarray(y, dtype=float))[0, 1])


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scenario, theta1), group in df.groupby(["scenario", "theta1"], sort=False):
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
                    "node_class": node_class,
                    "n_nodes": len(sub),
                    "pearson_r": pearsonr(sub["sim_centered"], sub["reconstructed_centered"]),
                    "rmse": float(np.sqrt(np.mean(residual**2))),
                    "mae": float(np.mean(np.abs(residual))),
                    "bias": float(np.mean(residual)),
                }
            )
    return pd.DataFrame(rows)


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.17,
        1.08,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        fontweight="bold",
    )


def set_equal_limits(ax, x: pd.Series, y: pd.Series) -> None:
    values = np.concatenate([np.asarray(x, dtype=float), np.asarray(y, dtype=float)])
    lo = float(np.nanpercentile(values, 0.5))
    hi = float(np.nanpercentile(values, 99.5))
    pad = 0.08 * max(hi - lo, 1.0)
    lo -= pad
    hi += pad
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    shared_ticks = ax.get_yticks()
    shared_ticks = shared_ticks[(shared_ticks >= lo) & (shared_ticks <= hi)]
    ax.set_xticks(shared_ticks)
    ax.set_yticks(shared_ticks)
    ax.plot([lo, hi], [lo, hi], color="#3b3b3b", lw=0.8, ls="--", zorder=1)


def plot_figure(df: pd.DataFrame, summary: pd.DataFrame, out_prefix: Path) -> None:
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

    fig, axes = plt.subplots(2, 3, figsize=(7.25, 4.85), constrained_layout=False)
    axes = axes.ravel()

    internal_color = "#b5b5b5"
    leaf_color = "#2b7c85"
    panels = list("ABCDEF")

    for i, (scenario, theta1) in enumerate(SCENARIOS):
        ax = axes[i]
        group = df[df["scenario"] == scenario].copy()
        internal = group[~group["is_leaf"].astype(bool)]
        leaves = group[group["is_leaf"].astype(bool)]

        ax.scatter(
            internal["sim_centered"],
            internal["reconstructed_centered"],
            s=9,
            color=internal_color,
            alpha=0.52,
            edgecolors="none",
            rasterized=True,
            label="Internal",
            zorder=2,
        )
        ax.scatter(
            leaves["sim_centered"],
            leaves["reconstructed_centered"],
            s=10,
            color=leaf_color,
            alpha=0.68,
            edgecolors="none",
            rasterized=True,
            label="Leaf",
            zorder=3,
        )
        set_equal_limits(ax, group["sim_centered"], group["reconstructed_centered"])
        ax.set_title(rf"OU ($\theta_1={theta1}$)")
        ax.set_xlabel("Simulated expression")
        ax.set_ylabel("Reconstructed expression")

        all_stats = summary[
            (summary["scenario"] == scenario) & (summary["node_class"] == "all nodes")
        ].iloc[0]
        ax.text(
            0.04,
            0.96,
            rf"$r$={all_stats['pearson_r']:.2f}" + "\n" + rf"MAE={all_stats['mae']:.2f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
        )
        add_panel_label(ax, panels[i])
        if i == 0:
            ax.legend(frameon=False, loc="lower right", handletextpad=0.2, borderpad=0.1)

    ax = axes[4]
    grouped = [df[df["scenario"] == scenario]["residual"].to_numpy() for scenario, _ in SCENARIOS]
    positions = np.arange(1, len(SCENARIOS) + 1)
    bp = ax.boxplot(
        grouped,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#202020", "linewidth": 1.0},
        boxprops={"linewidth": 0.8, "color": "#3b3b3b"},
        whiskerprops={"linewidth": 0.8, "color": "#3b3b3b"},
        capprops={"linewidth": 0.8, "color": "#3b3b3b"},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#d8e9ea")
    ax.set_xticks(positions)
    ax.set_xticklabels([rf"$\theta_1={theta1}$" for _, theta1 in SCENARIOS])
    ax.set_ylabel("Reconstructed - simulated")
    ax.set_title("Residuals")
    add_panel_label(ax, "E")

    ax = axes[5]
    all_summary = summary[summary["node_class"] == "all nodes"].copy()
    x = np.arange(len(all_summary))
    width = 0.36
    ax.bar(
        x - width / 2,
        all_summary["mae"],
        width=width,
        color="#2b7c85",
        label="MAE",
    )
    ax.bar(
        x + width / 2,
        all_summary["rmse"],
        width=width,
        color="#8d6e63",
        label="RMSE",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([rf"$\theta_1={theta1}$" for theta1 in all_summary["theta1"]])
    ax.set_ylabel("Expression error")
    ax.set_title("Accuracy summary")
    ax.legend(frameon=False, loc="upper right")
    add_panel_label(ax, "F")

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(length=2.5, width=0.7)

    fig.subplots_adjust(left=0.075, right=0.985, top=0.925, bottom=0.105, wspace=0.35, hspace=0.48)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".png"), dpi=450)
    plt.close(fig)


def write_caption(outdir: Path, out_prefix: Path, summary: pd.DataFrame, source_data: pd.DataFrame) -> None:
    all_stats = summary[summary["node_class"] == "all nodes"]
    stat_text = "; ".join(
        [
            rf"$\theta_1={int(row.theta1)}$: r={row.pearson_r:.2f}, MAE={row.mae:.2f}"
            for row in all_stats.itertuples(index=False)
        ]
    )
    source_files = source_data[["sim_history_file", "reconstruction_file"]].drop_duplicates()
    file_lines = "\n".join(
        [
            f"- `{row.sim_history_file}` and `{row.reconstruction_file}`"
            for row in source_files.itertuples(index=False)
        ]
    )
    caption = f"""# Supplementary figure: reconstructed versus simulated expression histories

Panels A-D compare LaVOUS posterior mean reconstructed expression with the simulated latent expression at matched tree nodes for the four Gene_500 OU examples shown in Fig. 2D. Values are centered on the simulated root for each scenario, matching the color centering used in Fig. 2D. Gray points are internal nodes and teal points are leaves; dashed lines indicate equality between simulated and reconstructed expression. Panel E shows residuals, defined as reconstructed minus simulated expression. Panel F summarizes mean absolute error (MAE) and root mean squared error (RMSE) across all nodes.

Summary across all nodes: {stat_text}.

Source files:
{file_lines}

Generated files:
- `{out_prefix.with_suffix('.pdf').name}`
- `{out_prefix.with_suffix('.png').name}`
- `source_data.tsv`
- `summary.tsv`
"""
    (outdir / "caption.md").write_text(caption)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--simdir", type=Path, default=Path("expression_simulation/simulation"))
    parser.add_argument("--reconstdir", type=Path, default=Path("expression_simulation/reconst"))
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("supplementary_figures/reconstruction_vs_simulated_expression"),
    )
    parser.add_argument("--prefix", default="lavous_reconstruction_vs_simulated_expression")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    source_data = pd.concat(
        [load_scenario(args.simdir, args.reconstdir, scenario, theta1) for scenario, theta1 in SCENARIOS],
        ignore_index=True,
    )
    summary = summarize(source_data)

    source_data.to_csv(args.outdir / "source_data.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)

    out_prefix = args.outdir / args.prefix
    plot_figure(source_data, summary, out_prefix)
    write_caption(args.outdir, out_prefix, summary, source_data)


if __name__ == "__main__":
    main()
