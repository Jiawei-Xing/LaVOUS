import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_curve


def load_lavous_scores(path):
    df = pd.read_csv(path, sep="\t")
    return df["lrt"].to_numpy(dtype=float)


def load_evogenex_scores(path):
    df = pd.read_csv(path)
    pvals = np.maximum(df["ou2_vs_ou1_pvalue"].iloc[::2].to_numpy(dtype=float), 1e-100)
    return -np.log10(pvals)


def draw_curve(ax, neg_scores, pos_scores, label, color):
    truth = np.array([False] * len(neg_scores) + [True] * len(pos_scores))
    scores = np.concatenate([neg_scores, pos_scores])
    fpr, tpr, _ = roc_curve(truth, scores)
    ax.plot(fpr, tpr, lw=1.5, color=color, label=f"{label} {auc(fpr, tpr):.3f}")


def plot_panel(ax, pos_label, neg_label):
    lavous_neg = load_lavous_scores(f"diff/diff_{neg_label}_chi-squared.tsv")
    lavous_pos = load_lavous_scores(f"diff/diff_{pos_label}_chi-squared.tsv")

    egx_neg_path = f"evogenex/egx_{neg_label}.csv"
    egx_pos_path = f"evogenex/egx_{pos_label}.csv"
    if os.path.exists(egx_neg_path) and os.path.exists(egx_pos_path):
        draw_curve(
            ax,
            load_evogenex_scores(egx_neg_path),
            load_evogenex_scores(egx_pos_path),
            "EvoGeneX",
            "tab:orange",
        )

    draw_curve(ax, lavous_neg, lavous_pos, "LaVOUS", "tab:green")

    ax.plot([0, 1], [0, 1], color="gray", lw=1, linestyle="--")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.legend(frameon=False, loc="lower right")


datasets = [
    [("Theta1", "Base"), ("t3r50", "rL"), ("t3r0", "r0")],
    [("t5r5", "Base"), ("t5r50", "rL"), ("t5r0", "r0")],
    [("t7r5", "Base"), ("t7r50", "rL"), ("t7r0", "r0")],
]

fig, axes = plt.subplots(3, 3, figsize=(9, 9))
for i, row in enumerate(datasets):
    for j, (pos, neg) in enumerate(row):
        plot_panel(axes[i, j], pos, neg)

fig.supxlabel("False Positive Rate")
fig.supylabel("True Positive Rate")
plt.tight_layout(pad=0.5)
fig.savefig("diff_ROC.png", dpi=300, bbox_inches="tight")
fig.savefig("diff_ROC.pdf", bbox_inches="tight")
plt.close(fig)
