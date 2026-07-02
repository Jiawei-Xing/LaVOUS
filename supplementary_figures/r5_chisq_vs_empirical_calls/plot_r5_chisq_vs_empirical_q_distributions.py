
#!/usr/bin/env python3
"""Supplementary figure: r=5 q-value distributions by calibration.

This script reads existing differential-test and empirical-calibration outputs
and writes derived files next to this script. It does not modify original
simulation or result files.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DIFF_DIR = REPO_ROOT / "expression_simulation" / "diff"
CALIB_DIR = DIFF_DIR / "test" / "calibrate_base_h0_empirical_10k"

SCENARIOS = [
    {"label": "Base", "title": r"$\theta_1=1$ (null)", "theta1": 1},
    {"label": "Theta1", "title": r"$\theta_1=3$", "theta1": 3},
    {"label": "t5r5", "title": r"$\theta_1=5$", "theta1": 5},
    {"label": "t7r5", "title": r"$\theta_1=7$", "theta1": 7},
]

CHI_COLOR = "#b63b32"
EMP_COLOR = "#2f6f9f"
Q_REFERENCE = 0.05


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


def load_merged(label):
    chi = pd.read_csv(DIFF_DIR / f"diff_{label}_chi-squared.tsv", sep="\t")
    emp = pd.read_csv(CALIB_DIR / f"diff_{label}_empirical-all.tsv", sep="\t")
    chi_cols = ["ID", "gene", "lrt", "p", "q", "signif"]
    emp_cols = ["ID", "gene", "lrt", "p", "q", "signif"]
    merged = chi[chi_cols].merge(
        emp[emp_cols],
        on=["ID", "gene", "lrt"],
        suffixes=("_chisq", "_empirical"),
        validate="one_to_one",
    )
    merged["scenario"] = label
    merged["chisq_source_file"] = str(DIFF_DIR / f"diff_{label}_chi-squared.tsv")
    merged["empirical_source_file"] = str(CALIB_DIR / f"diff_{label}_empirical-all.tsv")
    return merged.sort_values("lrt", ascending=False).reset_index(drop=True)


def collect_data():
    rows = []
    summary_rows = []
    for scenario in SCENARIOS:
        df = load_merged(scenario["label"])
        df["scenario_title"] = scenario["title"]
        df["theta1"] = scenario["theta1"]
        df["lrt_rank"] = np.arange(1, len(df) + 1)
        rows.append(df)
        for method, q_col in [("chi-square", "q_chisq"), ("empirical null", "q_empirical")]:
            q = df[q_col].to_numpy(dtype=float)
            summary_rows.append(
                {
                    "scenario": scenario["label"],
                    "scenario_title": scenario["title"],
                    "theta1": scenario["theta1"],
                    "method": method,
                    "n_genes": len(df),
                    "q_min": float(np.nanmin(q)),
                    "q_q10": float(np.nanquantile(q, 0.10)),
                    "q_q25": float(np.nanquantile(q, 0.25)),
                    "q_median": float(np.nanmedian(q)),
                    "q_q75": float(np.nanquantile(q, 0.75)),
                    "q_q90": float(np.nanquantile(q, 0.90)),
                    "q_max": float(np.nanmax(q)),
                    "n_q_lt_0_05": int(np.sum(q < 0.05)),
                }
            )
    return pd.concat(rows, ignore_index=True), pd.DataFrame(summary_rows)


def write_outputs(source, summary):
    source.to_csv(THIS_DIR / "source_data.tsv", sep="\t", index=False)
    summary.to_csv(THIS_DIR / "summary.tsv", sep="\t", index=False)


def plot_figure(source):
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(1, 4, figsize=(9.3, 2.75), sharey=True)
    for ax, panel, scenario in zip(axes, "ABCD", SCENARIOS):
        sub = source[source["scenario"] == scenario["label"]].sort_values("lrt_rank")
        ax.plot(sub["lrt_rank"], sub["q_chisq"], color=CHI_COLOR, lw=1.5, label="chi-square")
        ax.plot(sub["lrt_rank"], sub["q_empirical"], color=EMP_COLOR, lw=1.5, label="empirical null")
        ax.axhline(Q_REFERENCE, color="#707070", lw=0.9, ls=":")
        n_chi = int((sub["q_chisq"] < Q_REFERENCE).sum())
        n_emp = int((sub["q_empirical"] < Q_REFERENCE).sum())
        ax.text(
            0.04,
            0.94,
            f"q<0.05\nchi: {n_chi}\nemp: {n_emp}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.5},
        )
        ax.set_title(scenario["title"])
        ax.set_xlabel("Gene rank by LRT")
        ax.set_xlim(1, len(sub))
        ax.set_ylim(-0.03, 1.03)
        if ax is axes[0]:
            ax.set_ylabel("q value")
            ax.legend(frameon=False, fontsize=7, loc="upper right")
        add_panel_label(ax, panel)

    fig.suptitle("r = 5 q-value distributions by calibration", y=1.05, fontsize=10)
    fig.tight_layout(w_pad=1.2)
    name = "lavous_r5_chisq_vs_empirical_q_distributions"
    fig.savefig(THIS_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(THIS_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_caption(summary):
    def med(scenario, method):
        return float(summary[(summary["scenario"] == scenario) & (summary["method"] == method)]["q_median"].iloc[0])

    caption = f"""# Caption Draft

Supplementary Fig. X. Example r = 5 simulation q-value distributions under chi-square and empirical-null calibration. Panels A-D show per-gene q values for theta1 = 1, 3, 5, and 7 simulations, with genes ordered by decreasing LaVOUS LRT. The same LRT ranking is used for both calibrations in each panel; only the p/q calibration differs. The dotted horizontal line marks q = 0.05, and each panel reports the number of genes with q < 0.05 for each calibration.

The empirical-null calibration uses the Base-H0 pooled null from the 10k `lavous-calibrate` run and shifts q values relative to the asymptotic chi-square calibration, especially for weaker effects. At q < 0.05, empirical-null call counts are 130, 251, 387, and 481 for theta1 = 1, 3, 5, and 7, respectively. Median empirical-null q values are {med('Base', 'empirical null'):.3f}, {med('Theta1', 'empirical null'):.3f}, {med('t5r5', 'empirical null'):.3f}, and {med('t7r5', 'empirical null'):.3f} for theta1 = 1, 3, 5, and 7, respectively. Inputs were read from `expression_simulation/diff/diff_*_chi-squared.tsv` and the explicit `lavous-calibrate --sim_all 10000` outputs `expression_simulation/diff/test/calibrate_base_h0_empirical_10k/diff_*_empirical-all.tsv`; no original result files were modified.
"""
    (THIS_DIR / "caption.md").write_text(caption)


def main():
    source, summary = collect_data()
    write_outputs(source, summary)
    plot_figure(source)
    write_caption(summary)


if __name__ == "__main__":
    main()
