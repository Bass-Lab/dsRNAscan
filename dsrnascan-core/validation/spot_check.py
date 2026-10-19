"""Spot-check a whole-chromosome core run against dsRNAscan 0.5.x region runs.

A 0.5.x run with --start S (S - 1 a multiple of the step) has windows aligned with the whole-
chromosome run's, so in the region's interior (at least one window from either edge) the two must
report the same dsRNAs. Near the edges they may differ: the whole-chromosome run also has windows
that straddle the region's edges.

Usage: python spot_check.py <core whole-chromosome TSV> <0.5.x merged_results.txt> <start> <end> [margin]
"""
import sys

import pandas as pd

KEYS = ["Chromosome", "Strand", "i_start", "i_end", "j_start", "j_end"]
COLS = ["Score", "RawMatch", "Gaps", "dG(kcal/mol)", "base_pairs", "percent_paired", "longest_helix",
        "eff_i_start", "eff_i_end", "eff_j_start", "eff_j_end", "i_seq", "j_seq", "structure"]


def main():
    core_path, old_path, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    margin = int(sys.argv[5]) if len(sys.argv) > 5 else 10_000
    lo, hi = start + margin, end - margin
    inside = lambda d: d[(d.i_start >= lo) & (d.j_end <= hi)]
    core = inside(pd.read_csv(core_path, sep="\t"))
    old = inside(pd.read_csv(old_path, sep="\t"))
    m = old.merge(core, on=KEYS, how="outer", suffixes=("_o", "_c"), indicator=True)
    b = old.merge(core, on=KEYS, suffixes=("_o", "_c"))  # inner join: no NaN upcasting of integer columns
    diff = sum(int((b[c + "_o"].astype(str) != b[c + "_c"].astype(str)).sum()) for c in COLS)
    rows = int((b[[c + "_o" for c in COLS]].astype(str).values != b[[c + "_c" for c in COLS]].astype(str).values).any(axis=1).sum())
    print(f"interior {lo:,}-{hi:,}: 0.5.x region run {len(old)}, core whole chromosome {len(core)}, "
          f"identical coordinates {len(b)}, dsRNAs with any differing value {rows} ({diff} cells)")
    for _, r in m[m._merge != "both"].head(5).iterrows():
        print("   only in", "0.5.x" if r._merge == "left_only" else "core", r.Strand, r.i_start, r.i_end, r.j_start, r.j_end)


if __name__ == "__main__":
    main()
