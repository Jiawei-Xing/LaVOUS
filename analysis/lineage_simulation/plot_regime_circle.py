#!/usr/bin/env python
"""Plot a regime-colored tree using the simulation circular layout."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo


HERE = Path(__file__).resolve().parent
TREE_PATH = HERE / "tree.nwk"


def count_terminals(clade):
    if clade.is_terminal():
        return 1
    return sum(count_terminals(child) for child in clade.clades)


def ladderize(clade):
    if clade.is_terminal():
        return
    for child in clade.clades:
        ladderize(child)
    clade.clades.sort(key=count_terminals, reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=Path, default=TREE_PATH)
    parser.add_argument("--regime", type=Path, default=HERE / "regime.csv")
    parser.add_argument("--out-prefix", default="tree_regime")
    parser.add_argument("--title", default=None)
    parser.add_argument("--label-origins", action="store_true")
    args = parser.parse_args()

    out_png = HERE / f"{args.out_prefix}.png"
    out_pdf = HERE / f"{args.out_prefix}.pdf"

    tree = Phylo.read(args.tree, "newick")
    regimes = pd.read_csv(args.regime)
    regime_map = dict(zip(regimes["node_name"].astype(str), regimes["regime"].astype(str)))

    # Same circular layout as singlecellstochastics.stochas_sim.plot_sim_tree.
    ladderize(tree.root)
    all_nodes = list(tree.find_clades())
    leaves = list(tree.get_terminals())
    depths = tree.depths()
    angles = list(np.linspace(0, 2 * np.pi, len(leaves), endpoint=False))
    leaf_idx = [0]

    def assign_coords(node):
        node.r_coord = depths[node]
        if node.is_terminal():
            node.theta_coord = angles[leaf_idx[0]]
            leaf_idx[0] += 1
        else:
            for child in node.clades:
                assign_coords(child)
            node.theta_coord = np.mean([child.theta_coord for child in node.clades])

    assign_coords(tree.root)
    max_r = max(node.r_coord for node in all_nodes)

    parents = {}
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            parents[child] = parent

    origins = [
        child
        for child, parent in parents.items()
        if regime_map.get(child.name) == "1" and regime_map.get(parent.name) != "1"
    ]

    fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=200)
    ax.set_aspect("equal")
    ax.axis("off")

    tree_linewidth = 2.2 if len(leaves) > 500 else 3
    color_regime0 = "#9a9a9a"
    color_regime1 = "#d62728"

    def node_color(node):
        return color_regime1 if regime_map.get(node.name) == "1" else color_regime0

    def draw_tree(node):
        if not node.clades:
            return
        cr = node.r_coord
        col_parent = node_color(node)

        child_thetas = sorted([child.theta_coord for child in node.clades])
        if len(child_thetas) >= 2:
            t_min, t_max = child_thetas[0], child_thetas[-1]
            if t_max - t_min > np.pi:
                t_min, t_max = t_max, t_min + 2 * np.pi
            arc_t = np.linspace(t_min, t_max, max(20, int(abs(t_max - t_min) * 50)))
            ax.plot(
                cr * np.cos(arc_t),
                cr * np.sin(arc_t),
                color=col_parent,
                linewidth=tree_linewidth,
                solid_capstyle="round",
            )

        for child in node.clades:
            col_child = node_color(child)
            x0 = cr * np.cos(child.theta_coord)
            y0 = cr * np.sin(child.theta_coord)
            x1 = child.r_coord * np.cos(child.theta_coord)
            y1 = child.r_coord * np.sin(child.theta_coord)
            ax.plot(
                [x0, x1],
                [y0, y1],
                color=col_child,
                linewidth=tree_linewidth,
                solid_capstyle="butt",
            )
            draw_tree(child)

    draw_tree(tree.root)

    ax.scatter(
        [origin.r_coord * np.cos(origin.theta_coord) for origin in origins],
        [origin.r_coord * np.sin(origin.theta_coord) for origin in origins],
        color=color_regime1,
        s=28 if len(origins) <= 10 else 14,
        zorder=5,
    )

    # Label only independent regime origins when requested or when sparse enough
    # to keep the circular layout readable.
    if args.label_origins or len(origins) <= 8:
        for origin in sorted(origins, key=lambda node: node.theta_coord):
            label_r = max_r * 1.08
            lx = label_r * np.cos(origin.theta_coord)
            ly = label_r * np.sin(origin.theta_coord)
            ha = "left" if lx >= 0 else "right"
            ax.text(
                lx,
                ly,
                f"{origin.name} regime 1",
                color=color_regime1,
                fontsize=10,
                ha=ha,
                va="center",
                zorder=6,
            )

    n_regime1_tips = sum(regime_map.get(leaf.name) == "1" for leaf in leaves)
    title = args.title
    if title is None:
        title = (
            "Multi-regime tree: regime 1 in red "
            f"({n_regime1_tips}/{len(leaves)} tips, {len(origins)} origins)"
        )

    ax.set_title(title, fontsize=13, pad=12)
    margin = max_r * 1.22
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}")
    print(f"Wrote {out_pdf}")


if __name__ == "__main__":
    main()
