#!/usr/bin/env python3
"""
upgma.py

Supports two modes:
1) Sequence mode (default): hard‐coded DNA sequences → raw Hamming distances.
2) Matrix mode: read a user‐provided distance matrix (CSV or TSV with headers).

In both cases, runs UPGMA, prints distance matrices at each merge,
prints merge history + Newick string, and saves a dendrogram (with branch lengths)
to output/dendrogram.png.

Usage:
    # Default (sequence mode)
    python upgma.py

    # Matrix mode (provide a CSV/TSV file with first row=labels, first column=labels)
    python upgma.py --distances distances.csv
"""

import os
import sys
import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

# ─────────── 1. INPUT PARAMETERS ───────────
def parse_args():
    p = argparse.ArgumentParser(description="UPGMA clustering (from sequences or distance matrix)")
    p.add_argument(
        "--distances",
        "-d",
        metavar="FILE",
        help="Path to a CSV/TSV distance matrix (labels in first row/column)."
    )
    return p.parse_args()

# ─────────── 2. READ DISTANCE MATRIX ───────────
def read_distance_matrix(path):
    """
    Reads a distance matrix from a CSV/TSV file.
    Expect first row to be: , label1, label2, label3, ...
    Each subsequent row: label_i, d(i,label1), d(i,label2), ...
    Returns:
      - labels: list of strings
      - distances: dict mapping frozenset({i,j}) -> float distance
    """
    distances = {}
    with open(path, newline="") as f:
        # Sniff delimiter (comma or tab)
        sample = f.read(1024)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        reader = csv.reader(f, dialect)
        rows = list(reader)

    if not rows or len(rows) < 2:
        raise ValueError(f"Distance file '{path}' must have at least one header row and one data row.")

    headers = rows[0][1:]
    labels = headers.copy()

    # Build distances
    for row in rows[1:]:
        if len(row) < len(headers) + 1:
            raise ValueError("Each data row must have one leading label + same number of columns as header labels.")
        label1 = row[0]
        if label1 not in labels:
            raise ValueError(f"Label '{label1}' in first column not found among header labels.")
        for label2, val in zip(headers, row[1:]):
            try:
                d = float(val)
            except ValueError:
                raise ValueError(f"Cannot parse distance '{val}' as float for pair ({label1},{label2}).")
            # Only store once (i<j) via frozenset
            if label1 != label2:
                distances[frozenset({label1, label2})] = d
    return labels, distances

# ─────────── 3. DEFAULT SEQUENCES & DISTANCES ───────────
# Hard‐coded DNA sequences (example). Used if no --distances file provided.
sequences = {
    # exo 1
    "1" : "GTATAGGGGATATACTGAGAGCTATTTACA",
    "2" : "GTATTGGCGATATTCCGAGACCTATTTACT",
    "3" : "CTATTGGCCATATTCCGAGACCTATTTACT",
    "4" : "GTATAGGCGATATACCGAGACCTAATTACT",

    # exam 2020
    # "1" : "ACAAACAGTTCGATCGATTTGCAGTCTGGG",
    # "2" : "ACAAACAGTTTCTAGCGATTGCAGTCAGGG",
    # "3" : "ACAGACAGTTCGATCGATTTGCAGTCTCGG",
    # "4" : "ACTGACAGTTCGATCGATTTGCAGTCAGAG",
    # "5" : "ATTGACAGTTCGATCGATTTGCAGTCAGGA",
    # "6" : "TTTGACAGTTCGATCGATTTGCAGTCAGGG",

    # slides
    # "1" : "CATAGACCTGACGCCAGCTC",
    # "2" : "CATAGACCCGCCATGAGCTC",
    # "3" : "CGTAGACTGGGCGCCAGCTC",
    # "4" : "CCTAGACGTCGCGGCAGTCC",
}

def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Count the number of positions where seq1 and seq2 differ (raw count).
    """
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

def compute_initial_distances_from_sequences(seqs):
    """
    Build a dict of pairwise Hamming distances (raw counts) from given sequences.
    Returns:
      - labels: list of keys in seqs
      - distances: dict mapping frozenset({i,j}) -> float(raw_hamming_count)
    """
    labels = list(seqs.keys())
    distances = {}
    for i, j in combinations(labels, 2):
        raw = hamming_distance(seqs[i], seqs[j])
        distances[frozenset({i, j})] = float(raw)
    return labels, distances

# ─────────── 4. PRINT MATRIX ───────────
def print_matrix(distances, active_ids):
    """
    Nicely print the current “active” distance matrix.
    Rows/columns ordered by active_ids.
    """
    print("     " + "  ".join(f"{lab:>6}" for lab in active_ids))
    for i in active_ids:
        row = []
        for j in active_ids:
            if i == j:
                row.append("   0.0")
            else:
                d = distances.get(frozenset({i, j}), 0.0)
                # Print as integer if no fractional part, else one decimal place
                if d.is_integer():
                    row.append(f"{int(d):6d}")
                else:
                    row.append(f"{d:6.1f}")
        print(f"{i:>4} " + " ".join(row))

# ─────────── 5. UPGMA ALGORITHM ───────────
def upgma(labels, distances):
    """
    Perform UPGMA clustering given:
      - labels: list of leaf labels (strings)
      - distances: dict mapping frozenset({i,j}) -> float distance
    Returns:
      - clusters: dict cluster_id → {members, height, children}
      - root_id: the final cluster_id
      - merge_history: list of tuples (i, j, merge_distance)
    """
    # Initialize clusters: each label is its own cluster, height=0
    clusters = {lab: {"members": [lab], "height": 0.0, "children": []} for lab in labels}
    active_ids = labels.copy()
    merge_history = []

    print("\nInitial distance matrix:")
    print_matrix(distances, active_ids)

    next_cluster_id = 1 + max((int(x) for x in labels if x.isdigit()), default=0)
    next_cluster_id = str(next_cluster_id)

    while len(active_ids) > 1:
        # 1) Find closest pair
        pair_to_merge = min(distances, key=distances.get)
        i, j = tuple(pair_to_merge)
        min_dist = distances[pair_to_merge]
        new_id = str(next_cluster_id)
        new_height = min_dist / 2.0

        # 2) Create new cluster
        clusters[new_id] = {
            "members": clusters[i]["members"] + clusters[j]["members"],
            "height": new_height,
            "children": [
                (i, new_height - clusters[i]["height"]),
                (j, new_height - clusters[j]["height"]),
            ]
        }
        merge_history.append((i, j, min_dist))

        # 3) Update distances for new cluster
        for k in active_ids:
            if k not in (i, j):
                size_i = len(clusters[i]["members"])
                size_j = len(clusters[j]["members"])
                d_ik = distances.get(frozenset({i, k}), 0.0)
                d_jk = distances.get(frozenset({j, k}), 0.0)
                # Arithmetic mean of raw distances
                d_new_k = (size_i * d_ik + size_j * d_jk) / (size_i + size_j)
                distances[frozenset({new_id, k})] = d_new_k

        # 4) Remove old entries involving i or j
        to_remove = [key for key in distances if i in key or j in key]
        for key in to_remove:
            distances.pop(key)

        # 5) Update active_ids
        active_ids = [x for x in active_ids if x not in (i, j)]
        active_ids.append(new_id)
        next_cluster_id = str(int(next_cluster_id) + 1)

        # 6) Print updated matrix
        print(f"\nDistance matrix after merging {i} and {j} (distance = {min_dist:.1f}):")
        print_matrix(distances, active_ids)

    root = active_ids[0]
    return clusters, root, merge_history

# ─────────── 6. BUILD NEWICK ───────────
def build_newick(cluster_id, clusters):
    """
    Recursively build a Newick‐formatted string with branch lengths.
    """
    node = clusters[cluster_id]
    if not node["children"]:
        return f"{cluster_id}"
    parts = []
    for child_id, branch_len in node["children"]:
        subtree = build_newick(child_id, clusters)
        # format branch length as integer if whole, else one decimal
        if branch_len.is_integer():
            blabel = f"{int(branch_len)}"
        else:
            blabel = f"{branch_len:.1f}"
        parts.append(f"{subtree}:{blabel}")
    return f"({','.join(parts)})"

# ─────────── 7. DENDROGRAM COORDINATES ───────────
def get_leaves_order(cluster_id, clusters):
    """
    Return a list of leaf labels in left‐to‐right order for plotting.
    """
    node = clusters[cluster_id]
    if not node["children"]:
        return [cluster_id]
    left_id, _ = node["children"][0]
    right_id, _ = node["children"][1]
    return get_leaves_order(left_id, clusters) + get_leaves_order(right_id, clusters)

def assign_coords(cluster_id, clusters, leaf_positions, coords):
    """
    Recursively assign (x, y) to each node:
      - Leaves: y=0, x = leaf_positions[label]
      - Internals: x = midpoint of child x's, y = node height.
    """
    node = clusters[cluster_id]
    if not node["children"]:
        x = leaf_positions[cluster_id]
        y = 0.0
        coords[cluster_id] = (x, y)
        return x, y

    (l_id, _), (r_id, _) = node["children"]
    x_l, y_l = assign_coords(l_id, clusters, leaf_positions, coords)
    x_r, y_r = assign_coords(r_id, clusters, leaf_positions, coords)
    x_m = (x_l + x_r) / 2
    y_m = node["height"]
    coords[cluster_id] = (x_m, y_m)
    return x_m, y_m

# ─────────── 8. MAIN ───────────
if __name__ == "__main__":
    args = parse_args()

    # Ensure output folder
    os.makedirs("output", exist_ok=True)

    if args.distances:
        # Matrix mode
        try:
            labels, distances = read_distance_matrix(args.distances)
        except Exception as e:
            print(f"Error reading distance matrix: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Sequence mode: compute raw Hamming distances
        labels, distances = compute_initial_distances_from_sequences(sequences)

    # Run UPGMA
    clusters, root_id, history = upgma(labels, distances)

    # Print merge history
    print("\nMerge history (in order):")
    for idx, (a, b, d) in enumerate(history, 1):
        if d.is_integer():
            dist_label = f"{int(d)}"
            height_label = f"{int(d / 2)}"
        else:
            dist_label = f"{d:.1f}"
            height_label = f"{(d / 2):.1f}"
        print(f" {idx}. merged {a} and {b} at distance {dist_label} → height = {height_label}")

    # Print Newick string
    newick = build_newick(root_id, clusters) + ";"
    print("\nFinal tree in Newick format (with branch lengths):")
    print(newick)

    # Prepare for plotting
    leaves = get_leaves_order(root_id, clusters)
    leaf_positions = {leaf: idx for idx, leaf in enumerate(leaves)}
    coords = {}
    assign_coords(root_id, clusters, leaf_positions, coords)

    # Plot dendrogram + branch length labels
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title("UPGMA Dendrogram with Branch Lengths (Raw)")
    ax.set_ylabel("Distance (raw)")

    def plot_clade(cid):
        node = clusters[cid]
        x_c, y_c = coords[cid]
        if node["children"]:
            for child_id, br_len in node["children"]:
                x_ch, y_ch = coords[child_id]
                # Vertical line
                ax.plot([x_ch, x_ch], [y_ch, y_c], color="black")
                # Branch length label halfway up
                xm, ym = x_ch, (y_ch + y_c) / 2
                bl = int(br_len) if br_len.is_integer() else f"{br_len:.2f}"
                ax.text(xm, ym, f"{bl}", va="center", ha="right", fontsize=8)

            # Horizontal connector
            children_ids = [c[0] for c in node["children"]]
            xs = [coords[c][0] for c in children_ids]
            ax.plot([min(xs), max(xs)], [y_c, y_c], color="black")

            # Recurse
            for child_id, _ in node["children"]:
                plot_clade(child_id)

    plot_clade(root_id)

    # Label leaves along x-axis
    ax.set_xticks([leaf_positions[l] for l in leaves])
    ax.set_xticklabels(leaves)
    ax.invert_yaxis()
    plt.tight_layout()

    out_file = "output/dendrogram.png"
    plt.savefig(out_file, dpi=300)
    print(f"\nDendrogram image saved as: {out_file}")
