
#!/usr/bin/env python3
"""Supplementary figure: ELBO and LR quality by importance sampling.

The script reads existing fitted simulation outputs and writes derived figure
files next to this script. It does not modify original result files.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.special import gammaln
from scipy.stats import pearsonr, spearmanr

from singlecellstochastics.defaults import DEFAULT_LOG_S2_CLAMP
from singlecellstochastics.likelihood import ou_covariance_root_prior_torch
from singlecellstochastics.preprocess import process_data_OU
from singlecellstochastics.weights import theta_weight_W_torch


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[1]
DIFF_DIR = REPO_ROOT / "expression_simulation" / "diff"
DEFAULT_LABEL = "t5r5"
DEFAULT_N_GENES = 40
DEFAULT_N_SAMPLES = 512
DEFAULT_N_REPLICATES = 3
DEFAULT_BATCH_SIZE = 4
DEFAULT_PROPOSAL = "overdispersed"
DEFAULT_MIX = 0.8
DEFAULT_Q_SCALE = 2.0
DEFAULT_SEED = 20240625
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


def read_inputs(label, device):
    meta_path = DIFF_DIR / f"diff_{label}_meta.json"
    with meta_path.open() as handle:
        meta = json.load(handle)

    tree_files = [str(REPO_ROOT / meta["tree"])]
    expression_files = [str(REPO_ROOT / meta["expression"])]
    regime_files = [str(REPO_ROOT / meta["regime"])]
    library_files = [str(REPO_ROOT / meta["library"])] if meta.get("library") else [None]
    processed = process_data_OU(
        tree_files,
        expression_files,
        regime_files,
        library_files,
        meta["null"],
        device,
    )
    (
        _tree_list,
        cells_list,
        df_list,
        _diverge_list,
        _share_list,
        _epochs_list,
        _beta_list,
        diverge_list_torch,
        _share_list_torch,
        epochs_list_torch,
        beta_list_torch,
        regime_list,
        library_list,
    ) = processed

    diff = pd.read_csv(DIFF_DIR / f"diff_{label}_chi-squared.tsv", sep="\t")
    model_params = pd.read_csv(DIFF_DIR / f"diff_{label}_model-params.tsv", sep="\t")
    h0_q = pd.read_csv(DIFF_DIR / f"diff_{label}_h0_q-mean-std_0.tsv", sep="\t", index_col=0)
    h1_q = pd.read_csv(DIFF_DIR / f"diff_{label}_h1_q-mean-std_0.tsv", sep="\t", index_col=0)

    return {
        "meta": meta,
        "cells": cells_list[0],
        "counts": df_list[0],
        "epochs_list_torch": epochs_list_torch,
        "beta_list_torch": beta_list_torch,
        "diverge_list_torch": diverge_list_torch,
        "regimes": regime_list[0],
        "library": library_list[0],
        "diff": diff,
        "model_params": model_params,
        "h0_q": h0_q,
        "h1_q": h1_q,
    }


def select_genes(diff, n_genes):
    finite = diff[np.isfinite(diff["lrt"].to_numpy(dtype=float))].copy()
    finite = finite.sort_values("lrt").reset_index(drop=True)
    if n_genes >= len(finite):
        selected = finite.copy()
    else:
        idx = np.unique(np.round(np.linspace(0, len(finite) - 1, n_genes)).astype(int))
        selected = finite.iloc[idx].copy()
        if len(selected) < n_genes:
            missing = n_genes - len(selected)
            extra = finite.drop(selected.index).iloc[:missing]
            selected = pd.concat([selected, extra], axis=0).sort_values("lrt")
    selected["selection_rank"] = np.arange(1, len(selected) + 1)
    return selected


def q_tensor(q_df, cells, dtype, device):
    mean_cols = [f"q_mean_{cell}" for cell in cells]
    std_cols = [f"q_std_{cell}" for cell in cells]
    missing = [c for c in mean_cols + std_cols if c not in q_df.columns]
    if missing:
        raise ValueError(f"Missing q columns, first missing column: {missing[0]}")
    means = q_df.loc[:, mean_cols].to_numpy(dtype=np.float32)
    stds = q_df.loc[:, std_cols].to_numpy(dtype=np.float32)
    log_s2 = 2.0 * np.log(np.maximum(stds, 1e-12))
    q = np.concatenate([means, log_s2], axis=1)[:, None, :]
    return torch.tensor(q, dtype=dtype, device=device)


def fitted_params(model_params, genes, hypothesis, dtype, device):
    log_r_rows = []
    ou_rows = []
    for gene in genes:
        sub = model_params[(model_params["gene"] == gene) & (model_params["hypothesis"] == hypothesis)].copy()
        if sub.empty:
            raise ValueError(f"No {hypothesis} parameters found for {gene}")
        if hypothesis == "h0":
            row = sub.iloc[0]
            log_r_rows.append([math.log(float(row["r"]))])
            ou_rows.append([
                math.log(float(row["alpha"])),
                float(row["sigma"]),
                float(row["theta"]),
            ])
        elif hypothesis == "h1":
            sub["regime_str"] = sub["regime"].astype(str)
            sub = sub.sort_values("regime_str")
            first = sub.iloc[0]
            log_r_rows.append([math.log(float(first["r"]))])
            ou_rows.append(
                [math.log(float(first["alpha"])), float(first["sigma"])]
                + [float(x) for x in sub["theta"].to_numpy()]
            )
        else:
            raise ValueError(f"Unknown hypothesis: {hypothesis}")
    log_r = torch.tensor(log_r_rows, dtype=dtype, device=device)
    ou = torch.tensor(ou_rows, dtype=dtype, device=device)[:, None, :]
    return log_r, ou


def observation_constant(counts, lib_values):
    counts = np.asarray(counts, dtype=float)
    lib_values = np.asarray(lib_values, dtype=float)
    lib_values = np.maximum(lib_values, 1e-12)
    return np.sum(gammaln(counts + 1.0) - counts * np.log(lib_values[None, :]), axis=1)



def importance_log_weights(
    q_params,
    log_r,
    ou_params,
    mode,
    x_tensor,
    diverge_list_torch,
    epochs_list_torch,
    beta_list_torch,
    library_tensors,
    n_samples,
    proposal,
    mix,
    q_scale,
    nb=True,
    log_s2_clamp=DEFAULT_LOG_S2_CLAMP,
):
    """Return log importance weights for q, overdispersed q, or defensive prior mixtures."""
    if proposal not in {"q", "overdispersed", "defensive"}:
        raise ValueError(f"Unsupported proposal: {proposal}")
    if proposal == "q":
        mix = 1.0
    elif not 0.0 < mix < 1.0:
        raise ValueError(f"mix must be in (0, 1) for mixture proposals, got {mix}")
    if q_scale <= 1.0 and proposal == "overdispersed":
        raise ValueError(f"q_scale must be > 1 for overdispersed proposal, got {q_scale}")

    params = [q_params, log_r, ou_params]
    log_r_expanded = params[-2].unsqueeze(-1)
    alpha_log = params[-1][..., :1]
    other_ou = params[-1][..., 1:]
    log_weights_by_tree = []

    for tree_idx, x in enumerate(x_tensor):
        n_cells = x.shape[-1]
        q_mean = params[tree_idx][:, :, :n_cells]
        q_std = torch.exp(
            torch.clamp(
                params[tree_idx][:, :, n_cells : 2 * n_cells],
                log_s2_clamp[0],
                log_s2_clamp[1],
            )
            / 2.0
        )
        dist_q = torch.distributions.Normal(q_mean, q_std)
        samples_q = dist_q.sample((n_samples,))

        batch_size, n_sim, _n_cells = x.shape
        alpha = torch.exp(alpha_log[:, :, 0])
        sigma2 = other_ou[:, :, 0] ** 2
        thetas = other_ou[:, :, 1:]
        alpha_cov = alpha[:, :, None, None]
        sigma2_cov = sigma2[:, :, None, None]
        diverge = diverge_list_torch[tree_idx][None, None, :, :]
        cov = sigma2_cov * ou_covariance_root_prior_torch(alpha_cov, diverge)

        if mode == 1:
            mean = torch.ones(
                (batch_size, n_sim, n_cells), dtype=x.dtype, device=x.device
            ) * thetas[:, :, 0:1]
        elif mode == 2:
            W = theta_weight_W_torch(alpha, epochs_list_torch[tree_idx], beta_list_torch[tree_idx])
            mean = torch.matmul(W, thetas.unsqueeze(-1)).squeeze(-1)
        else:
            raise ValueError(f"Unsupported mode for OU IS: {mode}")

        dist_ou = torch.distributions.MultivariateNormal(loc=mean, covariance_matrix=cov)
        if proposal == "q":
            samples = samples_q
            log_p_z = dist_ou.log_prob(samples)
            denominator = dist_q.log_prob(samples).sum(dim=-1)
        elif proposal == "overdispersed":
            dist_wide = torch.distributions.Normal(q_mean, q_std * float(q_scale))
            samples_wide = dist_wide.sample((n_samples,))
            mask = (torch.rand(n_samples, batch_size, n_sim, device=x.device) < mix).to(x.dtype)
            samples = mask[..., None] * samples_q + (1.0 - mask[..., None]) * samples_wide
            log_p_z = dist_ou.log_prob(samples)
            log_q_z = dist_q.log_prob(samples).sum(dim=-1)
            log_wide_z = dist_wide.log_prob(samples).sum(dim=-1)
            denominator = torch.logsumexp(
                torch.stack([math.log(float(mix)) + log_q_z, math.log1p(-float(mix)) + log_wide_z], dim=0),
                dim=0,
            )
        else:
            samples_p = dist_ou.sample((n_samples,))
            mask = (torch.rand(n_samples, batch_size, n_sim, device=x.device) < mix).to(x.dtype)
            samples = mask[..., None] * samples_q + (1.0 - mask[..., None]) * samples_p
            log_p_z = dist_ou.log_prob(samples)
            log_q_z = dist_q.log_prob(samples).sum(dim=-1)
            denominator = torch.logsumexp(
                torch.stack([math.log(float(mix)) + log_q_z, math.log1p(-float(mix)) + log_p_z], dim=0),
                dim=0,
            )

        if nb:
            mu = torch.nn.functional.softplus(samples) * library_tensors[tree_idx]
            dist_count = torch.distributions.NegativeBinomial(
                total_count=torch.exp(log_r_expanded),
                logits=torch.log(mu) - log_r_expanded,
            )
        else:
            rate = torch.nn.functional.softplus(samples) * library_tensors[tree_idx]
            dist_count = torch.distributions.Poisson(rate=rate)
        log_p_x_given_z = dist_count.log_prob(x.unsqueeze(0)).sum(dim=-1)
        log_weights_by_tree.append(log_p_z + log_p_x_given_z - denominator)

    return torch.stack(log_weights_by_tree, dim=0).sum(dim=0)[:, :, 0]

def summarize_log_weights(log_weights):
    log_likelihood = torch.logsumexp(log_weights, dim=0) - math.log(float(log_weights.shape[0]))
    nll = -log_likelihood
    weights = torch.softmax(log_weights, dim=0)
    ess = 1.0 / torch.sum(weights**2, dim=0)

    shifted = torch.exp(log_weights - torch.max(log_weights, dim=0).values[None, :])
    mean_shifted = torch.mean(shifted, dim=0)
    sd_shifted = torch.std(shifted, dim=0, unbiased=True)
    loglik_se = (sd_shifted / math.sqrt(float(log_weights.shape[0]))) / torch.clamp(mean_shifted, min=1e-30)
    return nll, ess, loglik_se


def replicate_loglik_summary(rep_loglikes):
    if len(rep_loglikes) == 1:
        zeros = torch.zeros_like(rep_loglikes[0])
        nans = torch.full_like(rep_loglikes[0], float("nan"))
        return zeros, nans
    stacked = torch.stack(rep_loglikes, dim=0)
    sd = torch.std(stacked, dim=0, unbiased=True)
    se = sd / math.sqrt(float(len(rep_loglikes)))
    return se, sd


def compute_estimates(args):
    device = torch.device("cpu")
    dtype = torch.float32
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
    inputs = read_inputs(args.label, device)
    meta = inputs["meta"]
    if meta.get("root_mode", "stationary") != "stationary":
        raise ValueError("This figure script currently expects stationary OU root mode")

    selected = select_genes(inputs["diff"], args.n_genes)
    genes = selected["gene"].astype(str).tolist()
    counts_df = inputs["counts"]
    lib_values = inputs["library"].iloc[:, 0].to_numpy(dtype=float)
    constants = observation_constant(counts_df.loc[:, genes].to_numpy(dtype=float).T, lib_values)
    total_samples = int(args.n_samples * args.n_replicates)

    library_tensors = [torch.tensor(lib_values, dtype=dtype, device=device)]
    rows = []
    for start in range(0, len(genes), args.batch_size):
        end = min(start + args.batch_size, len(genes))
        batch_genes = genes[start:end]
        print(f"IS batch {start // args.batch_size + 1}: genes {start + 1}-{end} of {len(genes)}", flush=True)

        x_np = np.expand_dims(counts_df.loc[:, batch_genes].to_numpy(dtype=np.float32).T, axis=1)
        x_tensor = [torch.tensor(x_np, dtype=dtype, device=device)]

        h0_q = q_tensor(inputs["h0_q"].loc[batch_genes], inputs["cells"], dtype, device)
        h1_q = q_tensor(inputs["h1_q"].loc[batch_genes], inputs["cells"], dtype, device)
        h0_log_r, h0_ou = fitted_params(inputs["model_params"], batch_genes, "h0", dtype, device)
        h1_log_r, h1_ou = fitted_params(inputs["model_params"], batch_genes, "h1", dtype, device)

        h0_chunks = []
        h1_chunks = []
        h0_rep_loglikes = []
        h1_rep_loglikes = []
        with torch.no_grad():
            for rep in range(args.n_replicates):
                rep_seed = args.seed + 100000 * rep + start
                torch.manual_seed(rep_seed)
                h0_logw = importance_log_weights(
                    h0_q,
                    h0_log_r,
                    h0_ou,
                    1,
                    x_tensor,
                    inputs["diverge_list_torch"],
                    inputs["epochs_list_torch"],
                    inputs["beta_list_torch"],
                    library_tensors,
                    args.n_samples,
                    args.proposal,
                    args.mix,
                    args.q_scale,
                    nb=bool(meta.get("nb", True)),
                )
                torch.manual_seed(rep_seed)
                h1_logw = importance_log_weights(
                    h1_q,
                    h1_log_r,
                    h1_ou,
                    2,
                    x_tensor,
                    inputs["diverge_list_torch"],
                    inputs["epochs_list_torch"],
                    inputs["beta_list_torch"],
                    library_tensors,
                    args.n_samples,
                    args.proposal,
                    args.mix,
                    args.q_scale,
                    nb=bool(meta.get("nb", True)),
                )
                h0_chunks.append(h0_logw)
                h1_chunks.append(h1_logw)
                h0_rep_loglikes.append(torch.logsumexp(h0_logw, dim=0) - math.log(float(args.n_samples)))
                h1_rep_loglikes.append(torch.logsumexp(h1_logw, dim=0) - math.log(float(args.n_samples)))

            h0_all = torch.cat(h0_chunks, dim=0)
            h1_all = torch.cat(h1_chunks, dim=0)
            h0_nll, h0_ess, h0_loglik_se = summarize_log_weights(h0_all)
            h1_nll, h1_ess, h1_loglik_se = summarize_log_weights(h1_all)
            h0_rep_se, h0_rep_sd = replicate_loglik_summary(h0_rep_loglikes)
            h1_rep_se, h1_rep_sd = replicate_loglik_summary(h1_rep_loglikes)

        batch_diff = selected.set_index("gene").loc[batch_genes]
        for local_idx, gene in enumerate(batch_genes):
            row = batch_diff.loc[gene]
            const = float(constants[start + local_idx])
            h0_elbo = float(row["h0"])
            h1_elbo = float(row["h1"])
            h0_is_nll = float(h0_nll[local_idx])
            h1_is_nll = float(h1_nll[local_idx])
            is_lr = 2.0 * (h0_is_nll - h1_is_nll)
            elbo_lr = float(row["lrt"])
            min_ess = float(min(h0_ess[local_idx], h1_ess[local_idx]))
            rows.append(
                {
                    "label": args.label,
                    "n_samples_per_replicate": args.n_samples,
                    "n_replicates": args.n_replicates,
                    "total_importance_samples": total_samples,
                    "proposal": args.proposal,
                    "mix": args.mix,
                    "q_scale": args.q_scale,
                    "seed": args.seed,
                    "batch_size": args.batch_size,
                    "n_requested_genes": args.n_genes,
                    "selection_rank": int(row["selection_rank"]),
                    "ID": int(row["ID"]),
                    "gene": gene,
                    "elbo_h0_nll_no_constant": h0_elbo,
                    "elbo_h1_nll_no_constant": h1_elbo,
                    "observation_constant_added": const,
                    "elbo_h0_nll": h0_elbo + const,
                    "elbo_h1_nll": h1_elbo + const,
                    "is_h0_nll": h0_is_nll,
                    "is_h1_nll": h1_is_nll,
                    "elbo_lr": elbo_lr,
                    "is_lr": is_lr,
                    "is_minus_elbo_lr": is_lr - elbo_lr,
                    "h0_is_ess": float(h0_ess[local_idx]),
                    "h1_is_ess": float(h1_ess[local_idx]),
                    "min_is_ess": min_ess,
                    "min_is_ess_fraction": min_ess / float(total_samples),
                    "h0_is_loglik_mc_se": float(h0_loglik_se[local_idx]),
                    "h1_is_loglik_mc_se": float(h1_loglik_se[local_idx]),
                    "h0_replicate_loglik_se": float(h0_rep_se[local_idx]),
                    "h1_replicate_loglik_se": float(h1_rep_se[local_idx]),
                    "h0_replicate_loglik_sd": float(h0_rep_sd[local_idx]),
                    "h1_replicate_loglik_sd": float(h1_rep_sd[local_idx]),
                    "lr_replicate_mc_se": float(2.0 * torch.sqrt(h0_rep_se[local_idx] ** 2 + h1_rep_se[local_idx] ** 2)),
                    "p_chisq": float(row["p"]),
                    "q_chisq": float(row["q"]),
                    "signif_chisq": bool(row["signif"]),
                }
            )

    source = pd.DataFrame(rows).sort_values("selection_rank")
    source.to_csv(THIS_DIR / "source_data.tsv", sep="\t", index=False)
    source.to_csv(THIS_DIR / "importance_sampling_estimates.tsv", sep="\t", index=False)
    return source


def config_matches(source, args):
    if source.empty:
        return False
    checks = {
        "label": args.label,
        "n_samples_per_replicate": args.n_samples,
        "n_replicates": args.n_replicates,
        "total_importance_samples": args.n_samples * args.n_replicates,
        "proposal": args.proposal,
        "mix": float(args.mix),
        "q_scale": float(args.q_scale),
        "seed": args.seed,
        "batch_size": args.batch_size,
        "n_requested_genes": args.n_genes,
    }
    for column, expected in checks.items():
        if column not in source.columns:
            return False
        values = source[column].drop_duplicates()
        if len(values) != 1:
            return False
        observed = values.iloc[0]
        if isinstance(expected, float):
            if not np.isclose(float(observed), expected):
                return False
        elif observed != expected:
            return False
    return True


def load_or_compute(args):
    source_path = THIS_DIR / "source_data.tsv"
    if source_path.exists() and not args.force:
        source = pd.read_csv(source_path, sep="\t")
        if config_matches(source, args):
            print(f"Using cached IS estimates from {source_path}", flush=True)
            return source
    return compute_estimates(args)


def corr_metrics(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return pearsonr(x, y).statistic, spearmanr(x, y).statistic


def write_summary(df):
    h0_r, h0_rho = corr_metrics(df["elbo_h0_nll"], df["is_h0_nll"])
    h1_r, h1_rho = corr_metrics(df["elbo_h1_nll"], df["is_h1_nll"])
    lr_r, lr_rho = corr_metrics(df["elbo_lr"], df["is_lr"])
    summary = pd.DataFrame(
        [
            {
                "label": df["label"].iloc[0],
                "n_genes": len(df),
                "n_samples_per_replicate": int(df["n_samples_per_replicate"].iloc[0]),
                "n_replicates": int(df["n_replicates"].iloc[0]),
                "total_importance_samples": int(df["total_importance_samples"].iloc[0]),
                "proposal": df["proposal"].iloc[0],
                "mix": float(df["mix"].iloc[0]),
                "q_scale": float(df["q_scale"].iloc[0]),
                "pearson_elbo_vs_is_h0_nll": h0_r,
                "spearman_elbo_vs_is_h0_nll": h0_rho,
                "mae_h0_nll": np.mean(np.abs(df["is_h0_nll"] - df["elbo_h0_nll"])),
                "median_is_minus_elbo_h0_nll": np.median(df["is_h0_nll"] - df["elbo_h0_nll"]),
                "pearson_elbo_vs_is_h1_nll": h1_r,
                "spearman_elbo_vs_is_h1_nll": h1_rho,
                "mae_h1_nll": np.mean(np.abs(df["is_h1_nll"] - df["elbo_h1_nll"])),
                "median_is_minus_elbo_h1_nll": np.median(df["is_h1_nll"] - df["elbo_h1_nll"]),
                "pearson_elbo_lr_vs_is_lr": lr_r,
                "spearman_elbo_lr_vs_is_lr": lr_rho,
                "mae_lr": np.mean(np.abs(df["is_lr"] - df["elbo_lr"])),
                "median_is_minus_elbo_lr": np.median(df["is_lr"] - df["elbo_lr"]),
                "min_ess": np.min(df["min_is_ess"]),
                "median_min_ess": np.median(df["min_is_ess"]),
                "median_min_ess_fraction": np.median(df["min_is_ess_fraction"]),
                "median_lr_replicate_mc_se": np.nanmedian(df["lr_replicate_mc_se"]),
                "max_lr_replicate_mc_se": np.nanmax(df["lr_replicate_mc_se"]),
            }
        ]
    )
    summary.to_csv(THIS_DIR / "summary.tsv", sep="\t", index=False)
    return summary.iloc[0].to_dict()


def identity_limits(x, y, pad_fraction=0.04):
    values = np.concatenate([np.asarray(x, dtype=float), np.asarray(y, dtype=float)])
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    pad = max(1e-6, (hi - lo) * pad_fraction)
    return lo - pad, hi + pad


def lr_threshold_for_q(df, q_reference=Q_REFERENCE):
    label = str(df["label"].iloc[0])
    diff_path = DIFF_DIR / f"diff_{label}_chi-squared.tsv"
    if diff_path.exists():
        diff = pd.read_csv(diff_path, sep="\t")
        finite = diff[np.isfinite(diff["lrt"].to_numpy(dtype=float))]
        significant = finite[finite["q"] <= q_reference]
        if not significant.empty:
            return float(significant["lrt"].min())

    selected_significant = df[df["q_chisq"] <= q_reference]
    if not selected_significant.empty:
        return float(selected_significant["elbo_lr"].min())
    return None


def label_lr_threshold(ax, threshold, x_span, y_top):
    if threshold is None:
        return
    ax.axvline(threshold, color="#b63b32", lw=1.0, ls=":", zorder=1)
    ax.text(
        threshold + 0.015 * x_span,
        y_top,
        "q = 0.05\nELBO LR",
        rotation=90,
        color="#b63b32",
        ha="left",
        va="top",
        fontsize=7,
    )


def plot_figure(df, summary):
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(6.9, 5.45))
    axes = axes.ravel()
    lr_q_threshold = lr_threshold_for_q(df)
    signif_mask = df["q_chisq"] <= Q_REFERENCE

    panel_specs = [
        (axes[0], "A", "elbo_h0_nll", "is_h0_nll", "H0 marginal NLL", "#2f6f9f"),
        (axes[1], "B", "elbo_h1_nll", "is_h1_nll", "H1 marginal NLL", "#2ca25f"),
    ]
    for ax, panel, x_col, y_col, title, color in panel_specs:
        ax.scatter(df[x_col], df[y_col], s=22, color=color, alpha=0.78, linewidth=0)
        lo, hi = identity_limits(df[x_col], df[y_col])
        ax.plot([lo, hi], [lo, hi], color="#262626", lw=1, ls="--")
        r, rho = corr_metrics(df[x_col], df[y_col])
        mae = np.mean(np.abs(df[y_col] - df[x_col]))
        ax.text(
            0.04,
            0.96,
            f"r = {r:.2f}\nrho = {rho:.2f}\nMAE = {mae:.1f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=7,
        )
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_title(title)
        ax.set_xlabel("ELBO NLL")
        ax.set_ylabel("IS full NLL")
        add_panel_label(ax, panel)

    ax = axes[2]
    ax.scatter(
        df.loc[~signif_mask, "elbo_lr"],
        df.loc[~signif_mask, "is_lr"],
        s=23,
        color="#8f8f8f",
        alpha=0.78,
        linewidth=0,
        label="ELBO q > 0.05",
    )
    ax.scatter(
        df.loc[signif_mask, "elbo_lr"],
        df.loc[signif_mask, "is_lr"],
        s=23,
        color="#756bb1",
        alpha=0.82,
        linewidth=0,
        label="ELBO q <= 0.05",
    )
    _lo, hi = identity_limits(df["elbo_lr"], df["is_lr"], pad_fraction=0.07)
    lo = 0.0
    ax.plot([lo, hi], [lo, hi], color="#262626", lw=1, ls="--")
    if lr_q_threshold is not None:
        ax.axhline(lr_q_threshold, color="#b63b32", lw=0.9, ls=":", alpha=0.55, zorder=1)
        label_lr_threshold(ax, lr_q_threshold, hi - lo, lo + 0.52 * (hi - lo))
    ax.axhline(0, color="#bdbdbd", lw=0.8, zorder=0)
    ax.axvline(0, color="#bdbdbd", lw=0.8, zorder=0)
    ax.text(
        0.04,
        0.96,
        f"r = {summary['pearson_elbo_lr_vs_is_lr']:.2f}\nrho = {summary['spearman_elbo_lr_vs_is_lr']:.2f}\nMAE = {summary['mae_lr']:.1f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7,
    )
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_title("Likelihood-ratio statistic")
    ax.set_xlabel("ELBO LR")
    ax.set_ylabel("IS full LR")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    add_panel_label(ax, "C")


    ax = axes[3]
    lr_error = np.abs(df["is_lr"] - df["elbo_lr"]).to_numpy(dtype=float)
    sig_values = signif_mask.to_numpy()
    error_sets = [
        np.abs(df["is_h0_nll"] - df["elbo_h0_nll"]).to_numpy(dtype=float),
        np.abs(df["is_h1_nll"] - df["elbo_h1_nll"]).to_numpy(dtype=float),
        lr_error,
        lr_error[~sig_values],
    ]
    error_labels = ["H0 NLL", "H1 NLL", "LR", "LR\nq>0.05"]
    error_colors = ["#2f6f9f", "#2ca25f", "#756bb1", "#b63b32"]
    means = [float(np.mean(values)) for values in error_sets]
    x_pos = np.arange(len(error_sets))
    ax.bar(x_pos, means, color=error_colors, alpha=0.72, width=0.62, edgecolor="#262626", linewidth=0.7)
    for x_i, values in zip(x_pos, error_sets):
        offsets = np.linspace(-0.18, 0.18, len(values))
        ax.scatter(
            np.full(len(values), x_i) + offsets,
            np.sort(values),
            s=8,
            color="#262626",
            alpha=0.28,
            linewidth=0,
            zorder=3,
        )
    for x_i, value in zip(x_pos, means):
        ax.text(x_i, value + 0.55, f"{value:.1f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(error_labels)
    ax.set_ylim(0, max(max(np.max(v) for v in error_sets), max(means)) * 1.16)
    ax.set_title("Approximation error")
    ax.set_ylabel("Mean absolute error")
    add_panel_label(ax, "D")
    total_s = int(summary["total_importance_samples"])
    mix = float(summary["mix"])
    proposal = str(summary["proposal"])
    proposal_label = {"q": "q-proposal", "overdispersed": "overdispersed-q", "defensive": "defensive mixture"}.get(proposal, proposal)
    fig.suptitle(f"ELBO approximation quality by {proposal_label} IS likelihood (S={total_s}, mix={mix:g})", y=1.02, fontsize=10)
    fig.tight_layout(w_pad=1.5, h_pad=1.5)
    fig.savefig(THIS_DIR / "lavous_elbo_importance_sampling_quality.png", dpi=300, bbox_inches="tight")
    fig.savefig(THIS_DIR / "lavous_elbo_importance_sampling_quality.pdf", bbox_inches="tight")
    plt.close(fig)


def write_caption(df, summary):
    label = df["label"].iloc[0]
    n_genes = len(df)
    n_samples = int(summary["n_samples_per_replicate"])
    n_replicates = int(summary["n_replicates"])
    total_samples = int(summary["total_importance_samples"])
    mix = float(summary["mix"])
    proposal = str(summary["proposal"])
    q_scale = float(summary["q_scale"])
    noncalled = df[df["q_chisq"] > Q_REFERENCE]
    noncalled_lr_mae = float(np.mean(np.abs(noncalled["is_lr"] - noncalled["elbo_lr"])))
    n_noncalled = len(noncalled)
    if proposal == "q":
        proposal_sentence = "Importance sampling used the fitted variational distribution q(z|x) as the proposal"
        figure_phrase = "q-proposal importance-sampling"
    elif proposal == "overdispersed":
        proposal_sentence = f"Importance sampling used an overdispersed-q proposal, `{mix:.2f} q(z|x) + {1.0 - mix:.2f} q_{q_scale:g}(z|x)`, where q_{q_scale:g} has the same mean as q and {q_scale:g} times the posterior standard deviation"
        figure_phrase = "overdispersed-q importance-sampling"
    else:
        proposal_sentence = f"Importance sampling used a defensive proposal, `{mix:.2f} q(z|x) + {1.0 - mix:.2f} p(z)`"
        figure_phrase = "defensive importance-sampling"
    caption = f"""# Caption Draft

Supplementary Fig. X. ELBO approximation quality assessed by {figure_phrase} marginal likelihood estimates. Panels A and B compare the saved variational negative ELBOs with IS-estimated full marginal negative log likelihoods for H0 and H1, respectively, using a stratified subset of {n_genes} genes from the `{label}` simulation. Because the original differential-test run used `const=false`, the gene-wise negative-binomial observation constant was added back to the saved ELBO losses before these absolute NLL comparisons. Panel C compares the ELBO-based likelihood-ratio statistic, `2*(H0 - H1)`, with the corresponding IS-estimated full LR. The red dotted line marks the ELBO LR cutoff corresponding to q = 0.05 in the full `{label}` chi-square-calibrated differential-test result. Panel D compares mean absolute error for H0 NLL, H1 NLL, all LR values, and the LR subset with q > 0.05.

{proposal_sentence}, with {n_replicates} independent replicates of S = {n_samples} samples per gene and hypothesis ({total_samples} pooled samples total). The fitted variational distributions, OU parameters, and NB dispersions were read from `expression_simulation/diff/`; no model was re-optimized and no original result files were modified. Summary metrics for this subset: Pearson r = {summary['pearson_elbo_lr_vs_is_lr']:.3f} for ELBO LR versus IS full LR; mean absolute errors are H0 NLL = {summary['mae_h0_nll']:.3f}, H1 NLL = {summary['mae_h1_nll']:.3f}, LR = {summary['mae_lr']:.3f}, and LR among q > 0.05 genes = {noncalled_lr_mae:.3f} (n = {n_noncalled}); median replicate-based LR Monte Carlo SE = {summary['median_lr_replicate_mc_se']:.3f}.
"""
    (THIS_DIR / "caption.md").write_text(caption)


def parse_args():
    parser = argparse.ArgumentParser(description="Plot ELBO quality against IS-estimated full likelihood")
    parser.add_argument("--label", default=DEFAULT_LABEL, help="Simulation/diff label to use")
    parser.add_argument("--n-genes", type=int, default=DEFAULT_N_GENES, help="Number of genes stratified over ELBO LR")
    parser.add_argument("--n-samples", type=int, default=DEFAULT_N_SAMPLES, help="Importance samples per replicate, gene, and hypothesis")
    parser.add_argument("--n-replicates", type=int, default=DEFAULT_N_REPLICATES, help="Independent IS replicates to pool and use for MC stability")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="Genes per IS batch")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Torch RNG seed for IS samples")
    parser.add_argument("--proposal", choices=["q", "overdispersed", "defensive"], default=DEFAULT_PROPOSAL, help="Importance proposal family")
    parser.add_argument("--mix", type=float, default=DEFAULT_MIX, help="Weight on q(z|x) in mixture proposals")
    parser.add_argument("--q-scale", type=float, default=DEFAULT_Q_SCALE, help="Std multiplier for the overdispersed q component")
    parser.add_argument("--force", action="store_true", help="Recompute IS estimates even if cached source data match")
    return parser.parse_args()


def main():
    args = parse_args()
    df = load_or_compute(args)
    summary = write_summary(df)
    plot_figure(df, summary)
    write_caption(df, summary)
    print(f"Wrote {THIS_DIR / 'lavous_elbo_importance_sampling_quality.pdf'}", flush=True)


if __name__ == "__main__":
    main()
