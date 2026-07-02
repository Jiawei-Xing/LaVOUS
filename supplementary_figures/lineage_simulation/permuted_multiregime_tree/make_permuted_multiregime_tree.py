#!/usr/bin/env python
"""Build a strongly perturbed multi-origin regime lineage tree."""

from __future__ import annotations

import math
import random
from pathlib import Path

import pandas as pd
from Bio import Phylo


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
SOURCE_DIR = REPO_ROOT / "lineage_simulation"
TREE_PATH = SOURCE_DIR / "tree.nwk"

OUT_TREE = HERE / "tree_permuted_multiregime.nwk"
OUT_REGIME = HERE / "regime_permuted_multiregime.csv"
OUT_SWAPS = HERE / "permutation_summary.tsv"
OUT_REGIME_SUMMARY = HERE / "regime_summary.tsv"

RANDOM_SEED = 625

# Five disjoint origins from the existing multi-regime benchmark.  Together
# they keep the regime-1 tip fraction close to the original single-regime setup.
REGIME_ROOTS = ["c203", "c711", "c802", "c1054", "c1531"]

# Move each regime-origin clade onto a different background branch.  Regime
# labels are attached to clade objects before swapping, so the simulated regime
# history travels with the moved branch.
REGIME_SWAPS = [
    ("c203", "c447"),
    ("c711", "c852"),
    ("c802", "c1340"),
    ("c1054", "c1435"),
    ("c1531", "c1724"),
]

# Larger background-only moves make the non-regime portion visibly different.
LARGE_BACKGROUND_SWAPS = [
    ("c108", "c255"),
    ("c1273", "c819"),
    ("c324", "c271"),
]

N_BACKGROUND_SWAPS = 12


def parent_map(tree):
    parents = {}
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            parents[child] = parent
    return parents


def clade_by_name(tree):
    return {str(clade.name): clade for clade in tree.find_clades()}


def n_tips(clade):
    return len(clade.get_terminals())


def is_ancestor(ancestor, descendant):
    return descendant in set(ancestor.find_clades())


def all_regime_zero(clade):
    return all(getattr(desc, "regime", None) == "0" for desc in clade.find_clades())


def label_multiregime(tree):
    names = clade_by_name(tree)
    missing = [name for name in REGIME_ROOTS if name not in names]
    if missing:
        raise ValueError(f"Missing regime roots: {missing}")

    for clade in tree.find_clades():
        clade.regime = "0"
    for root_name in REGIME_ROOTS:
        for clade in names[root_name].find_clades():
            clade.regime = "1"


def swap_clades(tree, name_a, name_b, swap_type):
    names = clade_by_name(tree)
    clade_a = names[name_a]
    clade_b = names[name_b]
    parents = parent_map(tree)
    parent_a = parents[clade_a]
    parent_b = parents[clade_b]
    idx_a = parent_a.clades.index(clade_a)
    idx_b = parent_b.clades.index(clade_b)

    if clade_a is tree.root or clade_b is tree.root:
        raise ValueError("Cannot swap the root")
    if is_ancestor(clade_a, clade_b) or is_ancestor(clade_b, clade_a):
        raise ValueError(f"Cannot swap ancestor/descendant pair {name_a}, {name_b}")

    depths = tree.depths()
    depth_a = depths[clade_a]
    depth_b = depths[clade_b]
    new_branch_a = depth_a - depths[parent_b]
    new_branch_b = depth_b - depths[parent_a]
    if new_branch_a <= 0 or new_branch_b <= 0:
        raise ValueError(
            "Swap would require non-positive branch length: "
            f"{name_a}->{parent_b.name}={new_branch_a}, "
            f"{name_b}->{parent_a.name}={new_branch_b}"
        )

    row = {
        "swap_type": swap_type,
        "clade_a": name_a,
        "clade_b": name_b,
        "n_tips_a": n_tips(clade_a),
        "n_tips_b": n_tips(clade_b),
        "parent_a": str(parent_a.name),
        "parent_b": str(parent_b.name),
        "parent_depth_a": depths[parent_a],
        "parent_depth_b": depths[parent_b],
        "clade_depth_a": depth_a,
        "clade_depth_b": depth_b,
        "branch_length_a": clade_a.branch_length,
        "branch_length_b": clade_b.branch_length,
        "new_branch_length_a": new_branch_a,
        "new_branch_length_b": new_branch_b,
        "regime_a": clade_a.regime,
        "regime_b": clade_b.regime,
    }

    parent_a.clades[idx_a], parent_b.clades[idx_b] = clade_b, clade_a
    clade_a.branch_length = new_branch_a
    clade_b.branch_length = new_branch_b
    return row


def blocked_clades(tree):
    names = clade_by_name(tree)
    parents = parent_map(tree)
    blocked = set()
    for name_a, name_b in REGIME_SWAPS + LARGE_BACKGROUND_SWAPS:
        for name in (name_a, name_b):
            clade = names[name]
            blocked.update(clade.find_clades())
            parent = parents.get(clade)
            while parent is not None:
                blocked.add(parent)
                parent = parents.get(parent)
    return blocked


def choose_background_swaps(tree, rng):
    parents = parent_map(tree)
    depths = tree.depths()
    blocked = blocked_clades(tree)
    candidates = []
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root or clade.is_terminal() or clade in blocked:
            continue
        if not all_regime_zero(clade):
            continue
        tips = n_tips(clade)
        if tips < 8 or tips > 90:
            continue
        parent = parents[clade]
        candidates.append(
            {
                "name": str(clade.name),
                "clade": clade,
                "parent": parent,
                "parent_depth": depths[parent],
                "clade_depth": depths[clade],
                "n_tips": tips,
            }
        )

    pairs = []
    for i, item_a in enumerate(candidates):
        for item_b in candidates[i + 1 :]:
            clade_a = item_a["clade"]
            clade_b = item_b["clade"]
            if item_a["parent"] is item_b["parent"]:
                continue
            if is_ancestor(clade_a, clade_b) or is_ancestor(clade_b, clade_a):
                continue
            new_branch_a = item_a["clade_depth"] - item_b["parent_depth"]
            new_branch_b = item_b["clade_depth"] - item_a["parent_depth"]
            if new_branch_a <= 0 or new_branch_b <= 0:
                continue
            ratio = max(item_a["n_tips"], item_b["n_tips"]) / min(
                item_a["n_tips"], item_b["n_tips"]
            )
            if ratio > 2.5:
                continue
            parent_delta = abs(item_a["parent_depth"] - item_b["parent_depth"])
            if parent_delta > 4:
                continue
            size_score = abs(math.log(item_a["n_tips"] / item_b["n_tips"]))
            large_bonus = 0.04 * math.sqrt(min(item_a["n_tips"], item_b["n_tips"]))
            score = size_score + 0.12 * parent_delta - large_bonus
            pairs.append((score, item_a["name"], item_b["name"]))

    rng.shuffle(pairs)
    pairs.sort(key=lambda row: row[0])

    selected = []
    used = set()
    names = clade_by_name(tree)
    for _, name_a, name_b in pairs:
        clade_a = names[name_a]
        clade_b = names[name_b]
        if clade_a in used or clade_b in used:
            continue
        if any(is_ancestor(u, clade_a) or is_ancestor(clade_a, u) for u in used):
            continue
        if any(is_ancestor(u, clade_b) or is_ancestor(clade_b, u) for u in used):
            continue
        selected.append((name_a, name_b))
        used.add(clade_a)
        used.add(clade_b)
        if len(selected) == N_BACKGROUND_SWAPS:
            break

    return selected


def regime_summary(tree):
    parents = parent_map(tree)
    depths = tree.depths()
    rows = []
    regime1_tips = set()
    regime1_nodes = []
    regime1_branch_length = 0.0
    for clade in tree.find_clades(order="preorder"):
        if getattr(clade, "regime", None) != "1":
            continue
        regime1_nodes.append(clade)
        regime1_branch_length += float(clade.branch_length or 0.0)
        regime1_tips.update(tip.name for tip in clade.get_terminals())
        parent = parents.get(clade)
        if parent is None or getattr(parent, "regime", None) != "1":
            rows.append(
                {
                    "regime1_root": str(clade.name),
                    "n_tips": n_tips(clade),
                    "n_nodes": len(list(clade.find_clades())),
                    "depth": depths[clade],
                    "parent": str(parent.name) if parent else "",
                }
            )

    rows.append(
        {
            "regime1_root": "total",
            "n_tips": len(regime1_tips),
            "n_nodes": len(regime1_nodes),
            "depth": pd.NA,
            "parent": "",
            "regime1_branch_length": regime1_branch_length,
        }
    )
    return pd.DataFrame(rows)


def main():
    rng = random.Random(RANDOM_SEED)
    tree = Phylo.read(TREE_PATH, "newick")
    label_multiregime(tree)

    swap_rows = []
    for name_a, name_b in REGIME_SWAPS:
        swap_rows.append(swap_clades(tree, name_a, name_b, "regime"))

    for name_a, name_b in LARGE_BACKGROUND_SWAPS:
        swap_rows.append(swap_clades(tree, name_a, name_b, "background_large"))

    for name_a, name_b in choose_background_swaps(tree, rng):
        swap_rows.append(swap_clades(tree, name_a, name_b, "background"))

    Phylo.write(tree, OUT_TREE, "newick")
    pd.DataFrame(
        [
            {"node_name": str(clade.name), "regime": clade.regime}
            for clade in tree.find_clades(order="preorder")
        ]
    ).to_csv(OUT_REGIME, index=False)
    pd.DataFrame(swap_rows).to_csv(OUT_SWAPS, sep="\t", index=False)
    regime_summary(tree).to_csv(OUT_REGIME_SUMMARY, sep="\t", index=False)

    print(f"Wrote {OUT_TREE}")
    print(f"Wrote {OUT_REGIME}")
    print(f"Wrote {OUT_SWAPS}")
    print(f"Wrote {OUT_REGIME_SUMMARY}")


if __name__ == "__main__":
    main()
