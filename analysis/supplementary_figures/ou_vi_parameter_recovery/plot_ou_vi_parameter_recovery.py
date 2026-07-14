#!/usr/bin/env python3
"""Supplementary figure: OU and VI parameter recovery in simulations.

The script reads existing simulation and LaVOUS differential-test outputs and
writes all derived files next to this script. It does not modify original
results.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
SIM_DIR = REPO_ROOT / "expression_simulation" / "simulation"
DIFF_DIR = REPO_ROOT / "expression_simulation" / "diff"

SCENARIOS = [
    {"label": "Base", "title": "theta1=1", "theta0": 1.0, "theta1": 1.0, "alpha": 1.0, "sigma": 3.0, "color": "#3b6ea8"},
    {"label": "Theta1", "title": "theta1=3", "theta0": 1.0, "theta1": 3.0, "alpha": 1.0, "sigma": 3.0, "color": "#2ca25f"},
    {"label": "t5r5", "title": "theta1=5", "theta0": 1.0, "theta1": 5.0, "alpha": 1.0, "sigma": 3.0, "color": "#f28e2b"},
    {"label": "t7r5", "title": "theta1=7", "theta0": 1.0, "theta1": 7.0, "alpha": 1.0, "sigma": 3.0, "color": "#8f63b8"},
]
GENE_FOR_VI = "Gene_500"


def inv_softplus(y):
    y = np.asarray(y, dtype=float)
    y = np.maximum(y, 1e-12)
    return np.where(y > 20, y, np.log(np.expm1(y)))


def add_panel_label(ax, label):
    ax.text(
        -0.16,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="top",
    )


def load_model_params(label):
    path = DIFF_DIR / f"diff_{label}_model-params.tsv"
    return pd.read_csv(path, sep="\t")


def load_vi_params(label):
    path = DIFF_DIR / f"diff_{label}_h1_q-mean-std_0.tsv"
    return pd.read_csv(path, sep="\t", index_col=0)


def load_leaf_truth(label):
    path = SIM_DIR / f"sim_history_{label}.tsv"
    df = pd.read_csv(path, sep="\t")
    leaves = df[df["is_leaf"].astype(bool)].copy()
    # sim_history stores terminal expression after the simulator's softplus
    # transform, so invert it to compare with latent VI q means.
    leaves["true_latent"] = inv_softplus(leaves["sim_expr"].to_numpy(dtype=float))
    return leaves[["node_name", "sim_expr", "read_count", "true_latent"]]


def collect_ou_data():
    rows = []
    for scenario in SCENARIOS:
        label = scenario["label"]
        params = load_model_params(label)
        h1 = params[params["hypothesis"] == "h1"].copy()
        shared = h1.sort_values(["gene", "regime"]).drop_duplicates("gene")
        for _, row in shared.iterrows():
            inferred_stationary_var = float(row["sigma"]) ** 2 / (2.0 * float(row["alpha"]))
            true_stationary_var = float(scenario["sigma"]) ** 2 / (2.0 * float(scenario["alpha"]))
            rows.append(
                {
                    "scenario": label,
                    "scenario_title": scenario["title"],
                    "gene": row["gene"],
                    "parameter": "stationary_var",
                    "regime": "shared",
                    "inferred": inferred_stationary_var,
                    "truth": true_stationary_var,
                }
            )
        theta = h1.pivot(index="gene", columns="regime", values="theta")
        for gene, row in theta.iterrows():
            rows.append(
                {
                    "scenario": label,
                    "scenario_title": scenario["title"],
                    "gene": gene,
                    "parameter": "theta",
                    "regime": "0",
                    "inferred": float(row["0"]),
                    "truth": float(scenario["theta0"]),
                }
            )
            rows.append(
                {
                    "scenario": label,
                    "scenario_title": scenario["title"],
                    "gene": gene,
                    "parameter": "theta",
                    "regime": "1",
                    "inferred": float(row["1"]),
                    "truth": float(scenario["theta1"]),
                }
            )
    return pd.DataFrame(rows)


def collect_vi_data():
    rows = []
    for scenario in SCENARIOS:
        label = scenario["label"]
        q = load_vi_params(label)
        if GENE_FOR_VI not in q.index:
            raise ValueError(f"{GENE_FOR_VI} not found in {label} VI table")
        q_row = q.loc[GENE_FOR_VI]
        truth = load_leaf_truth(label).set_index("node_name")
        mean_cols = [c for c in q.columns if c.startswith("q_mean_")]
        for col in mean_cols:
            cell = col[len("q_mean_"):]
            std_col = "q_std_" + cell
            if cell not in truth.index or std_col not in q.columns:
                continue
            q_mean = float(q_row[col])
            q_std = float(q_row[std_col])
            true_latent = float(truth.loc[cell, "true_latent"])
            err = q_mean - true_latent
            rows.append(
                {
                    "scenario": label,
                    "scenario_title": scenario["title"],
                    "gene": GENE_FOR_VI,
                    "cell": cell,
                    "true_latent": true_latent,
                    "saved_terminal_expression": float(truth.loc[cell, "sim_expr"]),
                    "read_count": float(truth.loc[cell, "read_count"]),
                    "q_mean": q_mean,
                    "q_std": q_std,
                    "error": err,
                    "standardized_error": err / q_std if q_std > 0 else np.nan,
                    "covered_by_2sd": abs(err) <= 2.0 * q_std,
                }
            )
    return pd.DataFrame(rows)


def boxplot(ax, data, positions, color, widths=0.58):
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=widths,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.1},
        boxprops={"linewidth": 0.9},
        whiskerprops={"linewidth": 0.8},
        capprops={"linewidth": 0.8},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor(color)
        patch.set_alpha(0.72)
    return bp


def plot_figure(ou_df, vi_df):
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(8.1, 5.25))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.98, 1.08], hspace=0.43, wspace=0.36)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]
    titles = [s["title"] for s in SCENARIOS]
    tick_labels = [str(int(s["theta1"])) for s in SCENARIOS]
    x = np.arange(len(SCENARIOS))

    # A. Stationary variance recovery.
    ax = axes[0]
    data = [ou_df[(ou_df.scenario == s["label"]) & (ou_df.parameter == "stationary_var")]["inferred"].to_numpy() for s in SCENARIOS]
    boxplot(ax, data, x, "#9ecae1")
    true_var = SCENARIOS[0]["sigma"] ** 2 / (2.0 * SCENARIOS[0]["alpha"])
    ax.scatter(x, [true_var] * len(SCENARIOS), marker="D", s=20, color="#262626", zorder=4, label="truth")
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel(r"True $\theta_1$")
    ax.set_ylabel(r"Inferred $\sigma^2/(2\alpha)$")
    ax.set_title("OU stationary variance")
    add_panel_label(ax, "A")

    # B. Background optimum recovery.
    ax = axes[1]
    theta0 = [ou_df[(ou_df.scenario == s["label"]) & (ou_df.parameter == "theta") & (ou_df.regime == "0")]["inferred"].to_numpy() for s in SCENARIOS]
    boxplot(ax, theta0, x, "#bdbdbd")
    ax.scatter(x, [s["theta0"] for s in SCENARIOS], marker="D", s=20, color="#262626", zorder=4, label="truth")
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel(r"True $\theta_1$")
    ax.set_ylabel(r"Inferred $\theta_0$")
    ax.set_title("Background optimum")
    add_panel_label(ax, "B")

    # C. Selected-regime optimum recovery.
    ax = axes[2]
    theta1 = [ou_df[(ou_df.scenario == s["label"]) & (ou_df.parameter == "theta") & (ou_df.regime == "1")]["inferred"].to_numpy() for s in SCENARIOS]
    boxplot(ax, theta1, x, "#74c476")
    ax.scatter(x, [s["theta1"] for s in SCENARIOS], marker="D", s=20, color="#262626", zorder=4, label="truth")
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel(r"True $\theta_1$")
    ax.set_ylabel(r"Inferred $\theta_1$")
    ax.set_title("Selected-regime optimum")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    add_panel_label(ax, "C")

    # D. VI q mean recovery.
    ax = axes[3]
    min_v = min(vi_df["true_latent"].min(), vi_df["q_mean"].min())
    max_v = max(vi_df["true_latent"].max(), vi_df["q_mean"].max())
    pad = 0.04 * (max_v - min_v)
    for scenario in SCENARIOS:
        sub = vi_df[vi_df.scenario == scenario["label"]]
        ax.scatter(
            sub["true_latent"],
            sub["q_mean"],
            s=7,
            alpha=0.38,
            linewidth=0,
            color=scenario["color"],
            label=r"$\theta_1=%d$" % int(scenario["theta1"]),
            rasterized=True,
        )
    ax.plot([min_v - pad, max_v + pad], [min_v - pad, max_v + pad], color="#262626", linestyle="--", linewidth=1)
    corr = np.corrcoef(vi_df["true_latent"], vi_df["q_mean"])[0, 1]
    rmse = np.sqrt(np.mean(np.square(vi_df["error"])))
    ax.text(0.04, 0.96, f"r = {corr:.2f}\nRMSE = {rmse:.2f}", transform=ax.transAxes, va="top", ha="left", fontsize=7)
    ax.set_xlim(min_v - pad, max_v + pad)
    ax.set_ylim(min_v - pad, max_v + pad)
    ax.set_aspect("equal", adjustable="box")
    shared_ticks = ax.get_yticks()
    shared_ticks = shared_ticks[(shared_ticks >= min_v - pad) & (shared_ticks <= max_v + pad)]
    ax.set_xticks(shared_ticks)
    ax.set_yticks(shared_ticks)
    ax.set_xlabel(f"True latent leaf value ({GENE_FOR_VI})")
    ax.set_ylabel(r"VI $m_q$")
    ax.set_title("VI mean recovery")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right", markerscale=1.8)
    add_panel_label(ax, "D")

    # E. Standardized VI residuals.
    ax = axes[4]
    z = vi_df["standardized_error"].replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    bins = np.linspace(-5, 5, 51)
    ax.hist(z, bins=bins, density=True, color="#9ecae1", alpha=0.78, edgecolor="white", linewidth=0.25)
    xs = np.linspace(-5, 5, 300)
    normal = np.exp(-0.5 * xs**2) / np.sqrt(2 * np.pi)
    ax.plot(xs, normal, color="#262626", linestyle="--", linewidth=1.1, label="N(0,1)")
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_xlim(-5, 5)
    ax.set_xlabel(r"($m_q$ - truth) / $s_q$")
    ax.set_ylabel("Density")
    ax.set_title("VI standardized residuals")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    add_panel_label(ax, "E")

    # F. Two-standard-deviation coverage by scenario.
    ax = axes[5]
    coverage = vi_df.groupby("scenario")["covered_by_2sd"].mean().reindex([s["label"] for s in SCENARIOS])
    rmse_by = vi_df.groupby("scenario")["error"].apply(lambda e: np.sqrt(np.mean(np.square(e)))).reindex([s["label"] for s in SCENARIOS])
    colors = [s["color"] for s in SCENARIOS]
    ax.axhline(0.95, color="#b63b32", linestyle="--", linewidth=1.0, label="95%")
    for idx, (xi, val, err, color) in enumerate(zip(x, coverage.to_numpy(), rmse_by.to_numpy(), colors)):
        ax.scatter([xi], [val], s=44, color=color, edgecolor="#262626", linewidth=0.5, zorder=3)
        if idx == 0:
            ax.text(xi + 0.12, val - 0.006, f"{val:.2f}\nRMSE {err:.2f}", ha="left", va="top", fontsize=6.5)
        else:
            ax.text(xi, val + 0.004, f"{val:.2f}\nRMSE {err:.2f}", ha="center", va="bottom", fontsize=6.5)
    ax.set_ylim(0.90, 1.005)
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel(r"True $\theta_1$")
    ax.set_ylabel(r"Fraction within $m_q \pm 2s_q$")
    ax.set_title("VI uncertainty check")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    add_panel_label(ax, "F")

    fig.suptitle("Recovery of OU and variational parameters in r=5 OU simulations", y=1.02, fontsize=10)
    fig.savefig(THIS_DIR / "lavous_ou_vi_parameter_recovery.png", dpi=300, bbox_inches="tight")
    fig.savefig(THIS_DIR / "lavous_ou_vi_parameter_recovery.pdf", bbox_inches="tight")
    plt.close(fig)


def write_summary(ou_df, vi_df):
    ou_summary = (
        ou_df.groupby(["scenario", "scenario_title", "parameter", "regime", "truth"])["inferred"]
        .agg(n="count", median="median", q25=lambda x: np.quantile(x, 0.25), q75=lambda x: np.quantile(x, 0.75), mean="mean")
        .reset_index()
    )
    vi_summary = (
        vi_df.groupby(["scenario", "scenario_title"])
        .agg(
            n=("cell", "count"),
            corr=("q_mean", lambda x: np.corrcoef(vi_df.loc[x.index, "true_latent"], x)[0, 1]),
            rmse=("error", lambda x: np.sqrt(np.mean(np.square(x)))),
            bias=("error", "mean"),
            median_q_std=("q_std", "median"),
            coverage_2sd=("covered_by_2sd", "mean"),
        )
        .reset_index()
    )
    ou_summary.to_csv(THIS_DIR / "ou_parameter_summary.tsv", sep="\t", index=False)
    vi_summary.to_csv(THIS_DIR / "vi_parameter_summary.tsv", sep="\t", index=False)


def write_caption():
    caption = f"""# Caption Draft

Supplementary Fig. X. Recovery of OU and variational parameters in simulated data. Panels A-C summarize H1 OU parameter estimates across 500 genes from the r=5 OU differential-expression simulations (`Base`, `Theta1`, `t5r5`, and `t7r5`). Panel A shows the derived OU stationary variance, $\sigma^2/(2\alpha)$, rather than $\alpha$ and $\sigma$ separately. Black diamonds mark the simulation truth: $\sigma^2/(2\alpha)=4.5$, $\theta_0=1$, and $\theta_1 \in {{1, 3, 5, 7}}$. Panels D-F evaluate VI leaf-level parameters for `{GENE_FOR_VI}`, the gene for which `sim_history_*.tsv` stores the simulated expression history. Because the simulator writes terminal expression after the softplus transform, true terminal latent values were recovered with the inverse softplus before comparison with VI $m_q$. Panel E standardizes residuals by VI $s_q$, and Panel F shows the fraction of leaves covered by $m_q \pm 2s_q$.

Inputs were read from `expression_simulation/simulation/` and `expression_simulation/diff/`; no original result files were modified.
"""
    (THIS_DIR / "caption.md").write_text(caption)


def main():
    ou_df = collect_ou_data()
    vi_df = collect_vi_data()
    ou_df.to_csv(THIS_DIR / "ou_parameter_source_data.tsv", sep="\t", index=False)
    vi_df.to_csv(THIS_DIR / "vi_parameter_source_data.tsv", sep="\t", index=False)
    write_summary(ou_df, vi_df)
    write_caption()
    plot_figure(ou_df, vi_df)


if __name__ == "__main__":
    main()
