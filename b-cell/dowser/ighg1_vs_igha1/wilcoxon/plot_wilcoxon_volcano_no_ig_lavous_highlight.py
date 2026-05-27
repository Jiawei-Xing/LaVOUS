#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SCRIPT_DIR = Path(__file__).resolve().parent
DATA = SCRIPT_DIR.parent
WILCOXON = SCRIPT_DIR / "ighg1_vs_igha1_wilcoxon_1261_133.tsv"
LAVOUS = (
    DATA
    / "joint_outputs"
    / "diff"
    / "full"
    / "ighg1_vs_igha1_joint.full_empirical-each.tsv"
)
OUT_PREFIX = (
    SCRIPT_DIR / "volcano_ighg1_vs_igha1_wilcoxon_no_ig_lavous_highlight_1261_133"
)
IG_GENE_PATTERN = r"^(IGH|IGK|IGL)"
Q_THRESHOLD = 0.05
P_THRESHOLD = 0.05
X_CLIP = 8


LABEL_OFFSETS = {
    "XBP1": (-1.55, 0.9),
    "ACTG1": (0.85, 0.7),
    "FTL": (-1.7, 0.85),
    "RPS2": (0.85, 1.2),
}


def neg_log10(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    positive = values[values > 0]
    floor = positive.min() / 10 if not positive.empty else 1e-300
    return -np.log10(values.clip(lower=floor))


def main() -> None:
    wilcoxon = pd.read_csv(WILCOXON, sep="\t")
    lavous = pd.read_csv(LAVOUS, sep="\t")

    lavous = lavous.sort_values(["q", "p", "gene"], kind="mergesort").drop_duplicates(
        "gene", keep="first"
    )
    merged = wilcoxon.merge(
        lavous[["gene", "p", "q", "lrt", "signif"]].rename(
            columns={"p": "lavous_p", "q": "lavous_q", "signif": "lavous_signif"}
        ),
        left_on="symbol",
        right_on="gene",
        how="left",
        validate="many_to_one",
    )

    is_ig = merged["symbol"].astype(str).str.upper().str.match(IG_GENE_PATTERN)
    removed_ig = merged.loc[is_ig, "symbol"].drop_duplicates().sort_values().tolist()
    merged = merged.loc[~is_ig].copy()

    merged["wilcoxon_neg_log10_p"] = neg_log10(merged["pvalue"])
    merged["log2fc_clipped"] = merged["log2fc_norm_igha1_vs_ighg1"].clip(
        -X_CLIP, X_CLIP
    )
    merged["wilcoxon_bh_signif"] = merged["qvalue_bh"] < Q_THRESHOLD
    merged["lavous_q_signif"] = merged["lavous_q"] <= Q_THRESHOLD
    merged["lavous_direction"] = np.where(
        merged["log2fc_norm_igha1_vs_ighg1"] >= 0,
        "higher in IGHA1",
        "higher in IGHG1",
    )
    merged.to_csv(f"{OUT_PREFIX}.merged.tsv", sep="\t", index=False)

    background = merged[~merged["wilcoxon_bh_signif"] & ~merged["lavous_q_signif"]]
    wilcoxon_only = merged[merged["wilcoxon_bh_signif"] & ~merged["lavous_q_signif"]]
    lavous_sig = merged[merged["lavous_q_signif"]].copy()
    lavous_ighg1 = lavous_sig[lavous_sig["lavous_direction"] == "higher in IGHG1"]
    lavous_igha1 = lavous_sig[lavous_sig["lavous_direction"] == "higher in IGHA1"]

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
        }
    )
    fig, ax = plt.subplots(figsize=(10.8, 7.2), constrained_layout=True)

    ax.scatter(
        background["log2fc_clipped"],
        background["wilcoxon_neg_log10_p"],
        s=13,
        color="#d0d2d4",
        alpha=0.72,
        edgecolor="none",
        rasterized=True,
        zorder=1,
    )
    ax.scatter(
        wilcoxon_only["log2fc_clipped"],
        wilcoxon_only["wilcoxon_neg_log10_p"],
        s=30,
        color="#555861",
        alpha=0.95,
        edgecolor="white",
        linewidth=0.35,
        zorder=2,
    )
    ax.scatter(
        lavous_ighg1["log2fc_clipped"],
        lavous_ighg1["wilcoxon_neg_log10_p"],
        s=86,
        color="#0072b2",
        edgecolor="white",
        linewidth=0.8,
        zorder=4,
    )
    ax.scatter(
        lavous_igha1["log2fc_clipped"],
        lavous_igha1["wilcoxon_neg_log10_p"],
        s=86,
        color="#d55e00",
        edgecolor="white",
        linewidth=0.8,
        zorder=4,
    )

    ax.axvline(0, color="#8f8f8f", lw=1.1, ls="--", zorder=0)
    ax.axhline(-np.log10(P_THRESHOLD), color="#8da0b6", lw=1.0, ls=":", zorder=0)
    ax.text(
        X_CLIP - 0.05,
        -np.log10(P_THRESHOLD) + 0.04,
        "p = 0.05",
        ha="right",
        va="bottom",
        color="#667085",
        fontsize=10,
    )

    for _, row in lavous_sig.sort_values(
        "wilcoxon_neg_log10_p", ascending=False
    ).iterrows():
        dx, dy = LABEL_OFFSETS.get(row["symbol"], (0.8, 0.8))
        x = row["log2fc_clipped"]
        y = row["wilcoxon_neg_log10_p"]
        ha = "left" if dx >= 0 else "right"
        color = "#d55e00" if row["lavous_direction"] == "higher in IGHA1" else "#0072b2"
        symbol = row["symbol"]
        lavous_q = row["lavous_q"]
        wilcoxon_q = row["qvalue_bh"]
        label = f"{symbol}\nLAVOUS q={lavous_q:.3g}\nWilcoxon q={wilcoxon_q:.3g}"
        ax.annotate(
            label,
            xy=(x, y),
            xytext=(x + dx, y + dy),
            ha=ha,
            va="bottom",
            fontsize=9.5,
            color="#111827",
            bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": color, "lw": 0.9},
            arrowprops={"arrowstyle": "-", "color": color, "lw": 1.2},
            zorder=5,
        )

    max_y = max(merged["wilcoxon_neg_log10_p"].max(), -np.log10(P_THRESHOLD))
    ax.set_xlim(-X_CLIP - 0.35, X_CLIP + 0.35)
    ax.set_ylim(0, np.ceil(max_y + 0.8))
    ax.set_xlabel("log2 fold change, IGHA1 / IGHG1 (display clipped at +/-8)")
    ax.set_ylabel("-log10 Wilcoxon rank-sum p-value")
    ax.set_title("IGHG1 vs IGHA1 Wilcoxon differential expression (IG genes removed)")
    ax.text(
        0.01,
        0.985,
        "Negative: higher in IGHG1    Positive: higher in IGHA1",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color="#4b5563",
        fontsize=10.5,
    )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#d0d2d4",
            markeredgecolor="none",
            markersize=6,
            label=f"other non-IG genes (n={len(background)})",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#555861",
            markeredgecolor="white",
            markersize=7,
            label=f"Wilcoxon BH q < {Q_THRESHOLD:g} (n={len(wilcoxon_only)})",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#0072b2",
            markeredgecolor="white",
            markersize=8,
            label=f"LAVOUS q <= {Q_THRESHOLD:g}, IGHG1-high (n={len(lavous_ighg1)})",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#d55e00",
            markeredgecolor="white",
            markersize=8,
            label=f"LAVOUS q <= {Q_THRESHOLD:g}, IGHA1-high (n={len(lavous_igha1)})",
        ),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=10)

    fig.savefig(f"{OUT_PREFIX}.png", dpi=300)
    fig.savefig(f"{OUT_PREFIX}.svg")

    print(f"wrote\t{OUT_PREFIX}.png")
    print(f"wrote\t{OUT_PREFIX}.svg")
    print(f"wrote\t{OUT_PREFIX}.merged.tsv")
    print(f"non_ig_rows\t{len(merged)}")
    wilcoxon_bh_count = int(merged["wilcoxon_bh_signif"].sum())
    lavous_sig_count = int(merged["lavous_q_signif"].sum())
    print(f"wilcoxon_bh_q_lt_0.05\t{wilcoxon_bh_count}")
    print(f"lavous_q_le_0.05\t{lavous_sig_count}")
    print("lavous_sig\t" + ",".join(lavous_sig["symbol"].astype(str)))
    print(f"removed_ig_symbols\t{len(removed_ig)}")


if __name__ == "__main__":
    main()
