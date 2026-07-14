import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score


def load_scores(oup_file, egx_file):
    """Load raw ranking scores for each method (higher = more DE-like)."""
    df_oup = pd.read_csv(oup_file, sep="\t")
    df_egx = pd.read_csv(egx_file, sep=",")

    # OUP: LRT statistic
    score_oup = df_oup["lrt"].values

    # EGX: -log10(p), one row per gene (paired theta0/theta1 rows; p is on every row)
    p_egx = np.maximum(df_egx["ou2_vs_ou1_pvalue"].iloc[::2].values, 1e-100)
    score_egx = -np.log10(p_egx)

    return score_oup, score_egx


def plot_roc_curve(ax, oup_neg, oup_pos, egx_neg, egx_pos, title=None):
    s_oup_neg, s_egx_neg, = load_scores(oup_neg, egx_neg)
    s_oup_pos, s_egx_pos, = load_scores(oup_pos, egx_pos)

    n_neg = len(s_oup_neg)
    n_pos = len(s_oup_pos)
    truth = np.array([False] * n_neg + [True] * n_pos)

    for scores_neg, scores_pos, label, color in [
        (s_egx_neg, s_egx_pos, "EvoGeneX", "tab:orange"),
        (s_oup_neg, s_oup_pos, "LaVOUS", "tab:green")
    ]:
        p = np.concatenate((scores_neg, scores_pos))
        fpr, tpr, _ = roc_curve(truth, p)
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, lw=1.5, color=color, label=f"{label} {roc_auc:.3f}")

    ax.plot([0, 1], [0, 1], color="gray", lw=1, linestyle="--")
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    if title:
        ax.set_title(title)
    # Reorder legend: LaVOUS, EvoGeneX, SCOUT
    handles, labels = ax.get_legend_handles_labels()
    order = [0, 1]
    ax.legend([handles[i] for i in order], [labels[i] for i in order],
              frameon=False, loc="lower right")


# Define dataset pairs: (pos_label, neg_label)
datasets = [
    # Row 0
    [("Theta1", "Base"), ("t3r50", "rL"), ("t3r0", "r0")],
    # Row 1
    [("t5r5", "Base"), ("t5r50", "rL"), ("t5r0", "r0")],
    # Row 2
    [("t7r5", "Base"), ("t7r50", "rL"), ("t7r0", "r0")],
]

nrows = len(datasets)
ncols = len(datasets[0])
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))

for i, row in enumerate(datasets):
    for j, (pos, neg) in enumerate(row):
        oup_neg = f"diff/diff_{neg}_chi-squared.tsv"
        oup_pos = f"diff/diff_{pos}_chi-squared.tsv"
        egx_neg = f"evogenex/egx_{neg}.csv"
        egx_pos = f"evogenex/egx_{pos}.csv"
        plot_roc_curve(axes[i, j], oup_neg, oup_pos, egx_neg, egx_pos)

fig.supxlabel("False Positive Rate")
fig.supylabel("True Positive Rate")
plt.tight_layout(pad=0.5)
plt.savefig("diff_ROC.png", dpi=300, bbox_inches="tight")
plt.show()
