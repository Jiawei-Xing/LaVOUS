# ELBO Importance-Sampling Quality Figure

This folder is self-contained and read-only with respect to the original simulation outputs.
All generated files are written here unless a script variable is overridden.

## Regenerate the figure from existing fits

```bash
bash supplementary_figures/elbo_importance_sampling_quality/run_figure_importance_sampling.sh
```

Default settings use the `t5r5` simulation, 40 genes stratified over the saved ELBO LR,
3 independent IS replicates of 512 samples each, and an overdispersed-q proposal:
`0.8 q(z|x) + 0.2 q_2(z|x)`, where `q_2` has twice the fitted posterior standard deviation.

Useful overrides:

```bash
N_SAMPLES=1024 N_REPLICATES=3 bash supplementary_figures/elbo_importance_sampling_quality/run_figure_importance_sampling.sh
PROPOSAL=q MIX=1.0 bash supplementary_figures/elbo_importance_sampling_quality/run_figure_importance_sampling.sh
PROPOSAL=defensive MIX=0.95 bash supplementary_figures/elbo_importance_sampling_quality/run_figure_importance_sampling.sh
```

The defensive proposal mixes q with the OU prior and is much slower in this high-dimensional latent space.

## Run LAVOUS differential testing with package importance sampling

```bash
bash supplementary_figures/elbo_importance_sampling_quality/run_lavous_diff_with_importance_sampling.sh
```

This reruns `lavous-diff` into a separate output directory under this folder by default. The package CLI supports
q-proposal IS and q/prior defensive mixing through `--importance` and `--mix`; it does not implement the
overdispersed-q proposal used by the figure script.

Useful overrides:

```bash
OUTDIR=/tmp/lavous-t5r5-is IMPORTANCE=1024 MIX=1.0 BATCH=50 bash supplementary_figures/elbo_importance_sampling_quality/run_lavous_diff_with_importance_sampling.sh
```
