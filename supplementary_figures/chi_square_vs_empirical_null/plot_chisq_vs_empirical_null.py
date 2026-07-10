#!/usr/bin/env python3
"""Supplementary figure: null LRT QQ plots against chi-square.

This script is read-only with respect to existing simulation results. It writes
all derived figure files and source-data tables next to this script.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, kstest


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DIFF_DIR = REPO_ROOT / "expression_simulation" / "diff"
FDR_ALPHA = 0.05

# These are the three theta1 == theta0 null differential-test outputs present
# in expression_simulation/diff. Base is the r=5 baseline in run_diff.sh.
NULL_SETS = [
    ("Base", r"$\mathrm{r} = 5$", "#2f6f9f"),
    ("rL", r"$\mathrm{r} = 50$", "#2ca25f"),
    ("r0", r"$\mathrm{r} = \infty$", "#756bb1"),
]


def load_result(label):
    path = DIFF_DIR / f"diff_{label}_chi-squared.tsv"
    df = pd.read_csv(path, sep="\t")
    df["source_label"] = label
    return df


def finite_lrt(df):
    values = df["lrt"].to_numpy(dtype=float)
    return values[np.isfinite(values)]


def qq_data(values):
    observed = np.sort(np.asarray(values, dtype=float))
    n = len(observed)
    probs = (np.arange(1, n + 1) - 0.5) / n
    expected = chi2.ppf(probs, df=1)
    return expected, observed, probs


def add_panel_label(ax, label):
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
    )


def write_source_data(records):
    rows = []
    for label, title, _color, df, expected, observed, probs in records:
        sorted_df = df.loc[np.argsort(df["lrt"].to_numpy(dtype=float))].reset_index(drop=True)
        for rank, (_, row) in enumerate(sorted_df.iterrows(), start=1):
            rows.append(
                {
                    "source_label": label,
                    "panel_label": title,
                    "rank": rank,
                    "plotting_position": probs[rank - 1],
                    "chi_square_df1_quantile": expected[rank - 1],
                    "observed_lrt_quantile": observed[rank - 1],
                    "ID": row["ID"],
                    "gene": row["gene"],
                    "lrt": row["lrt"],
                    "p_chisq": row["p"],
                    "q_chisq": row["q"],
                    "signif_chisq": row["signif"],
                }
            )
    pd.DataFrame(rows).to_csv(THIS_DIR / "source_data.tsv", sep="\t", index=False)


def write_summary(records):
    rows = []
    chi95 = chi2.ppf(0.95, df=1)
    chi99 = chi2.ppf(0.99, df=1)
    for label, title, _color, df, _expected, _observed, _probs in records:
        lrt = finite_lrt(df)
        ks = kstest(lrt, chi2(df=1).cdf)
        rows.append(
            {
                "source_label": label,
                "panel_label": title,
                "n_genes": len(lrt),
                "lrt_min": np.min(lrt),
                "lrt_median": np.median(lrt),
                "lrt_q90": np.quantile(lrt, 0.90),
                "lrt_q95": np.quantile(lrt, 0.95),
                "lrt_q99": np.quantile(lrt, 0.99),
                "lrt_max": np.max(lrt),
                "chisq_df1_q95": chi95,
                "chisq_df1_q99": chi99,
                "fraction_lrt_ge_chisq_q95": np.mean(lrt >= chi95),
                "fraction_lrt_ge_chisq_q99": np.mean(lrt >= chi99),
                "p_le_0.05_calls": int((df["p"] <= FDR_ALPHA).sum()),
                "bh_q_le_0.05_calls": int((df["q"] <= FDR_ALPHA).sum()),
                "ks_statistic_vs_chisq_df1": ks.statistic,
                "ks_pvalue_vs_chisq_df1": ks.pvalue,
            }
        )
    pd.DataFrame(rows).to_csv(THIS_DIR / "summary.tsv", sep="\t", index=False)


def write_caption():
    caption = """# Caption Draft

Supplementary Fig. X. Null likelihood-ratio statistics compared with the asymptotic chi-square reference. QQ plots show the per-gene LaVOUS differential-test LRTs for the three theta1 = theta0 null simulation settings available in `expression_simulation/diff`: `Base` (r = 5 baseline), `rL` (r = 50), and `r0` (r = infinity). The black dashed line indicates equality with a chi-square distribution with one degree of freedom; red dotted vertical lines mark the expected-quantile position where genes first pass BH q = 0.05 in each panel. In all three null settings, the observed LRTs show strong upward tail deviation from the chi-square reference, explaining why chi-square-calibrated p-values are anti-conservative for these simulations.

No empirical null distributions were pooled across r settings; each panel uses only the corresponding null result file. Inputs were read from `expression_simulation/diff/`; no original result files were modified.
"""
    (THIS_DIR / "caption.md").write_text(caption)


def main():
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    records = []
    y_max = 0.0
    x_max = 0.0
    for label, title, color in NULL_SETS:
        df = load_result(label)
        lrt = finite_lrt(df)
        expected, observed, probs = qq_data(lrt)
        records.append((label, title, color, df, expected, observed, probs))
        y_max = max(y_max, float(np.max(observed)))
        x_max = max(x_max, float(np.max(expected)))

    fig, axes = plt.subplots(1, len(records), figsize=(7.6, 2.65), sharey=True)
    y_upper = y_max * 1.06
    y_lower = min(-1.0, min(float(np.min(record[5])) for record in records) * 1.08)
    diag_max = min(x_max * 1.04, y_upper)

    for ax, (panel, record) in zip(axes, zip("ABC", records)):
        label, title, color, df, expected, observed, _probs = record
        ax.scatter(expected, observed, s=10, alpha=0.68, linewidth=0, color=color)
        ax.plot([0, diag_max], [0, diag_max], color="#262626", lw=1, ls="--")
        q_sig = df[df["q"] <= FDR_ALPHA]
        if not q_sig.empty:
            q_threshold = float(q_sig["lrt"].min())
            boundary_idx = int(np.searchsorted(observed, q_threshold, side="left"))
            boundary_idx = min(max(boundary_idx, 0), len(expected) - 1)
            q_x = float(expected[boundary_idx])
            x_limit = x_max * 1.04
            label_x = q_x + 0.015 * x_limit
            label_ha = "left"
            if label_x > 0.96 * x_limit:
                label_x = q_x - 0.015 * x_limit
                label_ha = "right"
            ax.axvline(q_x, color="#b63b32", lw=0.9, ls=":")
            ax.text(
                label_x,
                y_upper * 0.92,
                "q = 0.05",
                rotation=90,
                color="#b63b32",
                ha=label_ha,
                va="top",
                fontsize=7,
            )
        ax.set_title(title)
        ax.set_xlim(0, x_max * 1.04)
        ax.set_ylim(y_lower, y_upper)
        ax.set_xlabel(r"Expected $\chi^2_1$ quantile")
        add_panel_label(ax, panel)

    axes[0].set_ylabel("Observed LRT quantile")
    fig.suptitle("Null differential-test LRTs vs chi-square reference", y=1.05, fontsize=10)
    fig.tight_layout(w_pad=1.0)
    fig.savefig(THIS_DIR / "lavous_chisq_vs_empirical_null.png", dpi=300, bbox_inches="tight")
    fig.savefig(THIS_DIR / "lavous_chisq_vs_empirical_null.pdf", bbox_inches="tight")
    plt.close(fig)

    write_source_data(records)
    write_summary(records)
    write_caption()


if __name__ == "__main__":
    main()
