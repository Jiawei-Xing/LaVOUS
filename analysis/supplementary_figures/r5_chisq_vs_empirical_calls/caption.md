# Caption Draft

Supplementary Fig. X. Example r = 5 simulation q-value distributions under chi-square and empirical-null calibration. Panels A-D show per-gene q values for theta1 = 1, 3, 5, and 7 simulations, with genes ordered by decreasing LaVOUS LRT. The same LRT ranking is used for both calibrations in each panel; only the p/q calibration differs. The dotted horizontal line marks q = 0.05.

The empirical-null calibration uses the Base-H0 pooled null from the 10k `lavous-calibrate` run and shifts q values relative to the asymptotic chi-square calibration, especially for weaker effects. Median empirical-null q values are 0.319, 0.049, 0.001, and 0.000 for theta1 = 1, 3, 5, and 7, respectively. Inputs were read from `expression_simulation/diff/diff_*_chi-squared.tsv` and the copied explicit `lavous-calibrate --sim_all 10000` outputs `supplementary_figures/r5_chisq_vs_empirical_calls/calibrate_base_h0_empirical_10k/diff_*_empirical-all.tsv`; no original result files were modified.
