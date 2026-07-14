"""
Nature-style 3×4 comparison figure.

Rows: Simulation | Reconstruction | Uncertainty
Cols: OU(θ₁=1) | OU(θ₁=3) | OU(θ₁=5) | OU(θ₁=7)

- Expression: diverging RdBu_r centered at 0 = each scenario's simulated root
  (white = root, blue = below, red = above). Sim and reconstruction rows are
  both re-centered on the *same* (true) root per column, so the colorbar reads
  as deviation from the true starting state.
- Variance: sequential Blues on log scale
- Shared colorbars

Usage:
    python plot_comparison.py --indir examples/reconst [--tree examples/input_data/tree.nwk]
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Bio import Phylo

# ── Config ───────────────────────────────────────────────────────
SCENARIOS = ["r0", "t3r0", "t5r0", "t7r0"]
COL_TITLES = [r"OU ($\theta_1$=1)", r"OU ($\theta_1$=3)",
              r"OU ($\theta_1$=5)", r"OU ($\theta_1$=7)"]
ROW_TITLES = ["Simulation", "Reconstruction", "Uncertainty"]
TITLE_FONTSIZE = 14
ROW_LABEL_FONTSIZE = 13
COLORBAR_LABEL_FONTSIZE = 12
COLORBAR_TICK_FONTSIZE = 10
# After re-centering, every scenario's root maps to 0 on the diverging cmap.
ROOT_VALUE = 0.0


# ── Tree layout ──────────────────────────────────────────────────

def layout_tree(newick_path):
    tree = Phylo.read(newick_path, "newick")

    # Ladderize (largest subtree first) so leaf order matches
    # stochas_sim.plot_sim_tree and reconstruct.plot_circular_tree.
    def count_terminals(clade):
        if clade.is_terminal():
            return 1
        return sum(count_terminals(c) for c in clade.clades)

    def ladderize(clade):
        if clade.is_terminal():
            return
        for child in clade.clades:
            ladderize(child)
        clade.clades.sort(key=count_terminals, reverse=True)

    ladderize(tree.root)

    all_nodes = list(tree.find_clades())
    leaves = list(tree.get_terminals())
    depths = tree.depths()
    angles = list(np.linspace(0, 2 * np.pi, len(leaves), endpoint=False))
    idx = [0]

    def assign(node):
        node.r_coord = depths[node]
        if node.is_terminal():
            node.theta_coord = angles[idx[0]]
            idx[0] += 1
        else:
            for c in node.clades:
                assign(c)
            node.theta_coord = np.mean([c.theta_coord for c in node.clades])

    assign(tree.root)
    max_r = max(n.r_coord for n in all_nodes)
    return tree, all_nodes, leaves, max_r


# ── Drawing functions ────────────────────────────────────────────

def draw_colored_tree(ax, tree, val_dict, cmap, norm, max_r,
                      outer_lens=None, outer_color_vals=None, branch_lw=2.0):
    """Draw tree with branches colored by val_dict.

    Optional outer ring bars at leaves: bar *length* scales with outer_lens
    (e.g. read counts), while bar *color* uses outer_color_vals on the
    same cmap/norm as branches. If omitted, bar color uses val_dict.
    """
    def draw(node):
        val = val_dict.get(node.name, np.nan)
        if np.isnan(val):
            val = ROOT_VALUE  # default to root for missing
        if node.clades:
            thetas = [c.theta_coord for c in node.clades]
            arc = np.linspace(min(thetas), max(thetas), 60)
            ax.plot(arc, [node.r_coord] * len(arc),
                    color=cmap(norm(val)), lw=branch_lw, zorder=1)
            for child in node.clades:
                c_val = val_dict.get(child.name, np.nan)
                if np.isnan(c_val):
                    c_val = val
                rs = np.linspace(node.r_coord, child.r_coord, 12)
                vs = np.linspace(val, c_val, 12)
                for k in range(len(rs) - 1):
                    ax.plot([child.theta_coord] * 2, [rs[k], rs[k + 1]],
                            color=cmap(norm(vs[k])),
                            lw=branch_lw, zorder=1, solid_capstyle='butt')
                draw(child)
    draw(tree.root)

    if outer_lens is not None:
        leaves = list(tree.get_terminals())
        len_vals = [outer_lens.get(l.name, 0) for l in leaves]
        len_max = max(len_vals) if len_vals else 1.0
        if len_max > 0:
            ax.fill_between(np.linspace(0, 2 * np.pi, 100),
                            max_r * 1.02, max_r * 1.22,
                            color='#f5f5f5', alpha=0.3)
            for leaf in leaves:
                ln = outer_lens.get(leaf.name, 0)
                if ln > 0:
                    bar_len = (ln / len_max) * (max_r * 0.20)
                    color_vals = outer_color_vals if outer_color_vals is not None else val_dict
                    leaf_val = color_vals.get(leaf.name, np.nan)
                    if np.isnan(leaf_val):
                        leaf_val = ROOT_VALUE
                    bar_color = cmap(norm(leaf_val))
                    ax.plot([leaf.theta_coord] * 2,
                            [max_r * 1.02, max_r * 1.02 + bar_len],
                            color=bar_color,
                            lw=2.5, solid_capstyle='butt')

    ax.axis('off')


# ── Data loading ─────────────────────────────────────────────────

def load_scenario(simdir, reconstdir, label):
    """Load all data for one scenario, re-centered on the simulated root.

    The first row of sim_history is the root (c0); we subtract its sim_expr
    from every node's sim_expr AND from every node's reconstructed mu so the
    sim and reconstruction rows share a common zero per column.
    """
    # Sim history (internal + leaf expression)
    sim_path = os.path.join(simdir, f"sim_history_{label}.tsv")
    if not os.path.exists(sim_path):
        print(f"  Warning: {sim_path} not found, simulation row will be blank")
        sim_dict = None
        sim_rc = None
        sim_root = 0.0  # nothing to subtract
    else:
        sdf = pd.read_csv(sim_path, sep='\t')
        sim_root = float(sdf['sim_expr'].iloc[0])  # root is row 0 (c0)
        centered = sdf['sim_expr'].astype(float) - sim_root
        sim_dict = dict(zip(sdf['node_name'].astype(str), centered))
        # read counts for outer ring (leaves only)
        leaf_rows = sdf[sdf['is_leaf']]
        sim_rc = dict(zip(leaf_rows['node_name'].astype(str),
                          leaf_rows['read_count']))

    # Read counts (leaves) — from readcounts file. Use the *last* gene
    # (Gene_<n_genes>), since that's the gene the reconstruction
    # (history_<label>_gene500.tsv) and sim_history were saved for. Using
    # the first column would show a different gene's counts in the
    # reconstruction row vs the simulation row.
    rc_path = os.path.join(simdir, f"readcounts_{label}.tsv")
    rc_df = pd.read_csv(rc_path, sep='\t', index_col=0)
    gene = rc_df.columns[-1]
    leaf_rc = dict(zip(rc_df.index.astype(str), rc_df[gene].astype(float)))

    # Reconstruction — re-center on the simulated root so sim and reconst
    # rows share a zero per column.
    rdf = pd.read_csv(os.path.join(reconstdir, f"history_{label}_gene500.tsv"),
                      sep='\t')
    mu_centered = rdf['infer_mu'].astype(float) - sim_root
    mu_dict = dict(zip(rdf['node_name'].astype(str), mu_centered))
    var_dict = dict(zip(rdf['node_name'].astype(str), rdf['infer_var']))

    return {
        'sim_dict': sim_dict,
        'sim_rc': sim_rc if sim_rc is not None else leaf_rc,
        'leaf_rc': leaf_rc,
        'mu_dict': mu_dict,
        'var_dict': var_dict,
        'sim_root': sim_root,
    }


# ── Main ─────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--simdir", default="expression_simulation/simulation")
    parser.add_argument("--reconstdir", default="expression_simulation/reconst")
    parser.add_argument("--tree", default="lineage_simulation/tree.nwk")
    args = parser.parse_args()

    # ── Load all scenarios ──
    data = {}
    all_sim, all_mu, all_var = [], [], []

    for label in SCENARIOS:
        try:
            d = load_scenario(args.simdir, args.reconstdir, label)
        except FileNotFoundError as e:
            print(f"Skipping {label}: {e}")
            continue
        data[label] = d
        if d['sim_dict'] is not None:
            all_sim.extend(d['sim_dict'].values())
        all_mu.extend(d['mu_dict'].values())
        all_var.extend(d['var_dict'].values())

    available = [s for s in SCENARIOS if s in data]
    n_cols = len(available)
    if n_cols == 0:
        print("No data found.")
        return

    # ── Color scales ──
    # Expression: RdBu_r centered at 0 (= simulated root, after re-centering).
    # Symmetric range from the 95th percentile of |deviation| across all
    # scenarios — ensures sim and reconstruction share one colorbar.
    all_dev = [v for v in all_sim + all_mu if np.isfinite(v)]
    expr_abs_max = np.percentile([abs(v) for v in all_dev], 95)
    expr_abs_max = max(expr_abs_max, 3.0)  # ensure minimum contrast
    norm_expr = mcolors.TwoSlopeNorm(vcenter=0.0,
                                      vmin=-expr_abs_max,
                                      vmax=expr_abs_max)
    cmap_expr = plt.cm.RdBu_r

    # Variance: log-scale Blues
    var_finite = [v for v in all_var if np.isfinite(v) and v > 0]
    norm_var = mcolors.LogNorm(vmin=min(var_finite), vmax=max(var_finite))
    cmap_var = plt.cm.Blues

    # ── Figure ──
    fig, axes = plt.subplots(3, n_cols,
                             subplot_kw={'projection': 'polar'},
                             figsize=(3.2 * n_cols, 3.2 * 3))
    if n_cols == 1:
        axes = axes[:, None]

    for j, label in enumerate(available):
        d = data[label]

        # Row 0: Simulation — colored tree + outer read-count bars
        tree, _, leaves, max_r = layout_tree(args.tree)
        if d['sim_dict'] is not None:
            draw_colored_tree(axes[0, j], tree, d['sim_dict'],
                              cmap_expr, norm_expr, max_r,
                              outer_lens=d['sim_rc'])
        else:
            axes[0, j].axis('off')

        # Row 1: Reconstruction — colored tree + outer read-count bars
        tree2, _, leaves2, max_r2 = layout_tree(args.tree)
        draw_colored_tree(axes[1, j], tree2, d['mu_dict'],
                          cmap_expr, norm_expr, max_r2,
                          outer_lens=d['leaf_rc'],
                          outer_color_vals=d['sim_dict'])

        # Row 2: Uncertainty — colored tree, no outer ring
        tree3, _, leaves3, max_r3 = layout_tree(args.tree)
        draw_colored_tree(axes[2, j], tree3, d['var_dict'],
                          cmap_var, norm_var, max_r3)

    # ── Column titles ──
    for j, label in enumerate(available):
        idx = SCENARIOS.index(label)
        axes[0, j].set_title(COL_TITLES[idx], fontsize=TITLE_FONTSIZE, fontweight='bold',
                             pad=12)

    # ── Row labels ──
    fig.canvas.draw()
    for i in range(3):
        bbox = axes[i, 0].get_position()
        fig.text(0.015, (bbox.y0 + bbox.y1) / 2, ROW_TITLES[i],
                 va='center', ha='center', fontsize=ROW_LABEL_FONTSIZE, fontweight='bold',
                 rotation=90)

    fig.subplots_adjust(left=0.05, right=0.92, bottom=0.06, top=0.94,
                        wspace=0.02, hspace=0.05)

    # ── Colorbars ──
    # Expression (rows 0-1)
    cax1 = fig.add_axes([0.94, 0.36, 0.013, 0.52])
    sm1 = plt.cm.ScalarMappable(norm=norm_expr, cmap=cmap_expr)
    cb1 = fig.colorbar(sm1, cax=cax1)
    cb1.set_label(r"Expression $-$ simulated root", fontsize=COLORBAR_LABEL_FONTSIZE)
    cb1.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)

    # Variance (row 2)
    cax2 = fig.add_axes([0.94, 0.06, 0.013, 0.22])
    sm2 = plt.cm.ScalarMappable(norm=norm_var, cmap=cmap_var)
    cb2 = fig.colorbar(sm2, cax=cax2)
    cb2.set_label("Variance", fontsize=COLORBAR_LABEL_FONTSIZE)
    cb2.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)

    out_path = os.path.join(args.reconstdir, "comparison.png")
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
