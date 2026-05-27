#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SCRIPT_DIR = Path(__file__).resolve().parent
DATA = SCRIPT_DIR.parent
LAVOUS = (
    DATA
    / "joint_outputs"
    / "diff"
    / "full"
    / "ighg1_vs_igha1_joint.full_empirical-each.tsv"
)
WILCOXON = SCRIPT_DIR / "ighg1_vs_igha1_wilcoxon_1261_133.tsv"
OUT_PREFIX = SCRIPT_DIR / "lavous_vs_wilcoxon_qvalues_no_ig_1261_133"
Q_THRESHOLD = 0.05
IG_GENE_PATTERN = r"^(IGH|IGK|IGL)"


def neg_log10_q(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    positive = values[values > 0]
    floor = positive.min() / 10 if not positive.empty else 1e-300
    return -np.log10(values.clip(lower=floor))


def quadrant(row: pd.Series) -> str:
    lav_sig = row["lavous_q"] <= Q_THRESHOLD
    wil_sig = row["wilcoxon_q"] <= Q_THRESHOLD
    if lav_sig and wil_sig:
        return "both"
    if lav_sig:
        return "lavous only"
    if wil_sig:
        return "wilcoxon only"
    return "neither"


def main() -> None:
    lavous = pd.read_csv(LAVOUS, sep="\t")
    wilcoxon = pd.read_csv(WILCOXON, sep="\t")
    wilcoxon = wilcoxon.sort_values(
        ["qvalue_bh", "pvalue", "symbol", "gene_id"], kind="mergesort"
    ).drop_duplicates("symbol", keep="first")

    merged = lavous[
        ["ID", "gene", "lrt", "p", "q", "h1_theta_IGHG1", "h1_theta_IGHA1"]
    ].merge(
        wilcoxon[
            [
                "symbol",
                "gene_id",
                "pvalue",
                "qvalue_bh",
                "diff_mean_log1p_igha1_minus_ighg1",
                "log2fc_norm_igha1_vs_ighg1",
            ]
        ],
        left_on="gene",
        right_on="symbol",
        how="left",
        validate="one_to_one",
    )
    missing = merged[merged["qvalue_bh"].isna()]
    if not missing.empty:
        genes = ", ".join(missing["gene"].astype(str))
        raise RuntimeError(f"Missing Wilcoxon q values for: {genes}")

    merged = merged.rename(
        columns={
            "p": "lavous_p",
            "q": "lavous_q",
            "pvalue": "wilcoxon_p",
            "qvalue_bh": "wilcoxon_q",
        }
    )
    is_ig = merged["gene"].astype(str).str.upper().str.match(IG_GENE_PATTERN)
    removed_ig = merged.loc[is_ig, "gene"].tolist()
    merged = merged.loc[~is_ig].copy()
    merged["lavous_neg_log10_q"] = neg_log10_q(merged["lavous_q"])
    merged["wilcoxon_neg_log10_q"] = neg_log10_q(merged["wilcoxon_q"])
    merged["category"] = merged.apply(quadrant, axis=1)
    merged["abs_delta_neg_log10_q"] = (
        merged["lavous_neg_log10_q"] - merged["wilcoxon_neg_log10_q"]
    ).abs()
    merged = merged.sort_values(
        ["lavous_q", "wilcoxon_q", "gene"], kind="mergesort"
    ).reset_index(drop=True)
    merged.to_csv(f"{OUT_PREFIX}.merged.tsv", sep="\t", index=False)

    colors = {
        "both": "#267c7c",
        "lavous only": "#c7502a",
        "wilcoxon only": "#6f55a4",
        "neither": "#7f7f7f",
    }
    order = ["neither", "wilcoxon only", "lavous only", "both"]

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
        }
    )
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)

    for cat in order:
        subset = merged[merged["category"] == cat]
        if subset.empty:
            continue
        ax.scatter(
            subset["lavous_neg_log10_q"],
            subset["wilcoxon_neg_log10_q"],
            s=52,
            c=colors[cat],
            edgecolor="white",
            linewidth=0.65,
            alpha=0.9,
            label=f"{cat} (n={len(subset)})",
            zorder=3,
        )

    max_axis = float(
        np.nanmax(
            [
                merged["lavous_neg_log10_q"].max(),
                merged["wilcoxon_neg_log10_q"].max(),
                -np.log10(Q_THRESHOLD),
            ]
        )
    )
    max_axis = np.ceil(max_axis + 0.5)
    ax.plot([0, max_axis], [0, max_axis], color="#333333", lw=1, ls=":", zorder=1)
    threshold = -np.log10(Q_THRESHOLD)
    ax.axvline(threshold, color="#555555", lw=1, ls="--", zorder=1)
    ax.axhline(threshold, color="#555555", lw=1, ls="--", zorder=1)

    label_genes = set(merged.nlargest(12, "abs_delta_neg_log10_q")["gene"].astype(str))
    label_genes.update(merged.nsmallest(5, "lavous_q")["gene"].astype(str))
    label_genes.update(merged.nsmallest(5, "wilcoxon_q")["gene"].astype(str))
    labels = merged[merged["gene"].astype(str).isin(label_genes)].copy()
    labels = labels.sort_values("wilcoxon_neg_log10_q", ascending=False).reset_index(
        drop=True
    )
    for i, row in labels.iterrows():
        dx = 0.08 if i % 2 == 0 else -0.08
        ha = "left" if dx > 0 else "right"
        ax.annotate(
            row["gene"],
            (row["lavous_neg_log10_q"], row["wilcoxon_neg_log10_q"]),
            xytext=(dx, 0.08),
            textcoords="offset fontsize",
            ha=ha,
            va="bottom",
            fontsize=8.5,
            color="#222222",
            arrowprops={"arrowstyle": "-", "lw": 0.35, "color": "#777777"},
        )

    ax.set_xlim(0, max_axis)
    ax.set_ylim(0, max_axis)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("LAVOUS empirical-each -log10(q)")
    ax.set_ylabel("Wilcoxon BH -log10(q)")
    ax.set_title("IGHG1 vs IGHA1 q-value comparison (IG genes removed)")
    ax.text(
        0.02,
        0.98,
        f"{len(merged)} non-IG LAVOUS-tested genes\nDashed lines: q = {Q_THRESHOLD}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="#333333",
    )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=colors[cat],
            markeredgecolor="white",
            markersize=7,
            label=f"{cat} (n={(merged['category'] == cat).sum()})",
        )
        for cat in order
        if (merged["category"] == cat).any()
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False)

    fig.savefig(f"{OUT_PREFIX}.png", dpi=300)
    fig.savefig(f"{OUT_PREFIX}.svg")
    print(f"wrote\t{OUT_PREFIX}.png")
    print(f"wrote\t{OUT_PREFIX}.svg")
    print(f"wrote\t{OUT_PREFIX}.merged.tsv")
    print(merged["category"].value_counts().to_csv(sep="\t", header=False).strip())
    print("removed_ig\t" + ",".join(removed_ig))


if __name__ == "__main__":
    main()
