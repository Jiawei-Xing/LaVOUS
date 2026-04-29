# LAVOUS Source Map

This document maps the publication workflows to the code paths in
`singlecellstochastics/`. It is intended to make review and maintenance easier
without changing the mathematical implementation.

## Core Data Flow

1. Command-line entry points parse files and optimization options.
2. `preprocess.py` reads trees, count matrices, library-size factors, and regime
   labels. It normalizes branch lengths, aligns cells to tree leaves, and builds
   covariance helper matrices.
3. `lrt.py`, `plasticity.py`, or `selection.py` initialize variational and model
   parameters for each gene batch.
4. `optimize.py` optimizes the ELBO or direct Gaussian likelihood.
5. `elbo.py`, `likelihood.py`, `approx.py`, and `weights.py` evaluate the model.
6. `output.py` writes test statistics and FDR-adjusted results.

## Workflow Entry Points

### Expression Heritability

CLI:

- `lavous-heritability`
- `run-plasticity-test`

Main function:

- `plasticity.run_plasticity_test`

Model comparison:

- H0: BM/NB with Pagel lambda fixed at 0.
- H1: BM/NB with free Pagel lambda.

Important code:

- `preprocess.process_data_BM`
- `optimize.Lq_optimize_torch_BM`
- `likelihood.bm_neg_log_lik_torch_kkt`

### Differential Expression

CLI:

- `lavous-diff`
- `run-diff-test`

Main functions:

- `ou_diff.run_diff_test`
- `lrt.likelihood_ratio_test`

Model comparison:

- H0: OU/NB with one shared theta.
- H1: OU/NB with regime-specific theta values.

Important code:

- `preprocess.process_data_OU`
- `optimize.Lq_optimize_torch_OU`
- `elbo.Lq_neg_log_lik_torch`
- `likelihood.ou_neg_log_lik_torch_kkt`
- `output.save_result`
- `output.output_results`
- `output.output_model_params`

### Empirical-Null Calibration

CLI:

- `lavous-calibrate`
- `run-calibrate`

Main function:

- `calibrate.run_calibrate`

Purpose:

- Reuses H0 parameters and metadata from a completed differential-expression
  run.
- Simulates null read-count datasets.
- Refits the LRT on simulated data and replaces chi-squared p-values with
  empirical p-values.

Important code:

- `simulate.simulate_null_all`
- `simulate.simulate_null_each`
- `lrt.likelihood_ratio_test`

### History Reconstruction

CLI:

- `lavous-reconstruct`
- `run-reconst`

Main function:

- `reconstruct.run_reconst`

Purpose:

- Uses fitted OU/BM transition parameters and leaf variational beliefs to run
  Gaussian belief propagation from leaves to root and back down the tree.

Important code:

- `reconstruct.upward_pass`
- `reconstruct.downward_pass`
- `reconstruct.plot_circular_tree`

## Parameter Conventions

Variational parameters are stored internally as:

- first half: latent means
- second half: `log(sigma_q^2)`

Use `qparam.export_mean_std_tensor` before writing files so outputs remain
mean-plus-standard-deviation tables.

Differential-expression result files report:

- `delta_nll`: `h0 - h1`, the optimized negative-log-likelihood difference.
- `lrt`: `2 * delta_nll`, the chi-squared likelihood-ratio statistic.

Empirical-null calibration compares `delta_nll` values because the factor of two
does not change the tail probability.

OU model tensors use log alpha internally and report `alpha = exp(log_alpha)`.
Negative-binomial dispersion uses log r internally and reports `r = exp(log_r)`.

Pagel lambda modes:

- `0`: fixed lambda = 0, star-tree covariance.
- `1`: fixed lambda = 1, original tree covariance.
- `2`: optimized lambda constrained to `(0, 1)` by sigmoid.

## Generated and Non-Source Files

The following paths are ignored by `.gitignore` and should not be considered
publication source:

- `__pycache__/`
- `.ipynb_checkpoints/`
- `wandb/`
- test and simulation output directories matching `test*`
- build and egg-info directories
