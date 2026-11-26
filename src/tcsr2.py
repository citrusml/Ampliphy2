#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tcsr2.py
========
Taxonomic Congruence Score without Root v2 (TCSr2 / TCSr_new) calculator.

Implements the unrooted graph-based clade splitting penalty and normalized
similarity congruence score:
    TCSr2 = 1.0 - (TCSr_raw / max_penalty)
in [0.0, 1.0], where 1.0 indicates perfect taxonomic congruence.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Union

import numpy as np
import toytree


def build_unrooted_graph(
    tree_newick: str,
) -> Tuple[Dict[int, Set[int]], Dict[str, int], int]:
    """
    Parse a Newick tree and construct an unrooted undirected graph adjacency map.

    Parameters
    ----------
    tree_newick : str
        Newick tree string. (unitless)

    Returns
    -------
    adj : Dict[int, Set[int]]
        Undirected adjacency list representing the tree graph. (unitless)
    leaf_map : Dict[str, int]
        Mapping from clean leaf name to node ID. (unitless)
    root_id : int
        Original root node ID. (unitless)
    """
    t = toytree.tree(tree_newick.strip())
    adj: Dict[int, Set[int]] = {}
    node_id_map: Dict[Any, int] = {}
    leaf_map: Dict[str, int] = {}

    for idx, node in enumerate(t.treenode.traverse()):
        node_id_map[node] = idx
        adj[idx] = set()
        if node.is_leaf():
            clean_name = str(node.name).strip("'\"")
            leaf_map[clean_name] = idx

    for node in t.treenode.traverse():
        u = node_id_map[node]
        for child in node.children:
            v = node_id_map[child]
            adj[u].add(v)
            adj[v].add(u)

    # Collapse artificial degree-2 root node to preserve true unrooted topology
    root_id = node_id_map[t.treenode]
    if len(adj.get(root_id, set())) == 2:
        nbrs = list(adj[root_id])
        adj[nbrs[0]].remove(root_id)
        adj[nbrs[1]].remove(root_id)
        adj[nbrs[0]].add(nbrs[1])
        adj[nbrs[1]].add(nbrs[0])
        del adj[root_id]

    return adj, leaf_map, root_id


def count_clade_groups_in_tree(
    adj: Dict[int, Set[int]],
    leaf_ids: List[int],
    clade_leaf_ids: List[int],
) -> int:
    """
    Compute the number of disjoint groups (connected components) of clade leaves
    in the unrooted tree graph.

    Parameters
    ----------
    adj : Dict[int, Set[int]]
        Undirected adjacency list of the tree graph. (unitless)
    leaf_ids : List[int]
        Node IDs of all leaves in the tree. (unitless)
    clade_leaf_ids : List[int]
        Node IDs of leaves possessing the evaluated clade. (unitless)

    Returns
    -------
    int
        Number of groups k_c formed by the clade on the tree. (count)
    """
    n_leaves = len(leaf_ids)
    n_c = len(clade_leaf_ids)

    if n_c == n_leaves:
        return 0
    if n_c <= 1 or (n_leaves - n_c) <= 1:
        return 1

    non_c_leaves = set(leaf_ids) - set(clade_leaf_ids)
    degrees = {u: len(nbrs) for u, nbrs in adj.items()}
    pruned: Set[int] = set(clade_leaf_ids)
    queue = list(clade_leaf_ids)

    idx = 0
    while idx < len(queue):
        u = queue[idx]
        idx += 1
        for v in adj[u]:
            if v not in pruned:
                degrees[v] -= 1
                if v not in non_c_leaves and degrees[v] <= 1:
                    pruned.add(v)
                    queue.append(v)

    boundary_count = 0
    for u in pruned:
        for v in adj[u]:
            if v not in pruned:
                boundary_count += 1

    return boundary_count


def compute_tree_tcsr2(
    tree_newick: str,
    taxa_map: Dict[str, Set[str]],
    species_code_len: int = 5,
) -> float:
    """
    Compute TCSr2 = 1.0 - (TCSr_raw / max_penalty) for an unrooted Newick tree.

    Parameters
    ----------
    tree_newick : str
        Input Newick tree string. (unitless)
    taxa_map : Dict[str, Set[str]]
        Mapping of species code prefix (e.g., 'HUMAN') to set of lineage clades. (unitless)
    species_code_len : int, default=5
        Length of mnemonic prefix in leaf name. (count)

    Returns
    -------
    tcsr2 : float
        Unrooted taxonomic congruence score in [0.0, 1.0], or NaN on invalid input. (unitless)
    """
    if not tree_newick or not tree_newick.strip():
        return float("nan")

    try:
        adj, leaf_map, _ = build_unrooted_graph(tree_newick)
    except Exception:
        return float("nan")

    num_leaves = len(leaf_map)
    if num_leaves < 3:
        return 1.0

    leaf_lineages: Dict[int, Set[str]] = {}
    all_clades: Set[str] = set()
    for leaf_name, leaf_id in leaf_map.items():
        code = leaf_name[:species_code_len]
        lin = taxa_map.get(code, set())
        leaf_lineages[leaf_id] = lin
        all_clades.update(lin)

    clade_to_leaves: Dict[str, List[int]] = {c: [] for c in all_clades}
    for leaf_id, lin in leaf_lineages.items():
        for c in lin:
            clade_to_leaves[c].append(leaf_id)

    leaf_ids = list(leaf_map.values())
    raw_score = 0
    max_penalty = 0

    for _c, c_leaf_ids in clade_to_leaves.items():
        n_c = len(c_leaf_ids)
        if n_c == num_leaves or n_c <= 1 or (num_leaves - n_c) <= 1:
            continue
        k_c = count_clade_groups_in_tree(adj, leaf_ids, c_leaf_ids)
        penalty = max(0, k_c - 1)
        raw_score += penalty
        max_penalty += (n_c - 1)

    return 1.0 - (float(raw_score) / max_penalty) if max_penalty > 0 else 1.0


def load_taxa_lineages(lineage_path: Union[str, Path]) -> Dict[str, Set[str]]:
    """
    Load species taxonomic lineage mapping from an AmpliPhy 2-column lineage.tsv.

    Parameters
    ----------
    lineage_path : Union[str, Path]
        Path to lineage.tsv (query<TAB>lineage). (unitless)

    Returns
    -------
    taxa_map : Dict[str, Set[str]]
        Dictionary mapping 5-letter species code to set of taxonomic clades. (unitless)
    """
    taxa_map: Dict[str, Set[str]] = {}
    with open(lineage_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("query"):
                continue
            parts = s.split("\t")
            if len(parts) >= 2:
                code = parts[0].strip()[:5]
                taxa_map[code] = {x.strip() for x in parts[1].split(",") if x.strip()}
    return taxa_map


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute TCSr2 for a Newick tree")
    parser.add_argument("tree", help="Input Newick tree file")
    parser.add_argument("lineage", help="Lineage TSV file")
    args = parser.parse_args()

    tree_str = Path(args.tree).read_text(encoding="utf-8").strip()
    taxa_map = load_taxa_lineages(args.lineage)
    score = compute_tree_tcsr2(tree_str, taxa_map)
    if np.isnan(score):
        print("nan")
    else:
        print(f"{score:.6f}")


if __name__ == "__main__":
    main()
