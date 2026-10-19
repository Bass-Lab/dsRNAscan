"""End-to-end comparison with the dsRNAscan 0.5.x pipeline on the same sequence.

Old: windows + einverted + RNAduplex + filters + nested elimination, run beforehand with
     python dsrnascan/dsRNAscan.py <fa> --no-ml --output-dir <dir> -c 8
New: dsrnascan-core end to end (search, Rust duplexfold, same filters, box rule), in each of its
     modes: default (0.5.x-equivalent settings), --compact, --banded.

Usage: python vs_pipeline.py <fa> <old *_merged_results.txt>
"""
import csv
import sys
from collections import Counter

from common import ov, run_core

MODES = {"default (0.5.x replay)": [], "--banded": ["--banded"],
         "--windowless": ["--windowless"], "--windowless --compact": ["--windowless", "--compact"]}


def recip(a, b, side):
    s, e = ("i_start", "i_end") if side == 5 else ("j_start", "j_end")
    o = ov(a[s], a[e], b[s], b[e])
    return o >= 0.5 * (a[e] - a[s] + 1) and o >= 0.5 * (b[e] - b[s] + 1)


def same_dsrna(a, b):
    return a["strand"] == b["strand"] and recip(a, b, 5) and recip(a, b, 3)


def contains(big, small):
    """Both arms of `small` are >= 90% inside the matching arms of `big`, on the same strand."""
    return (big["strand"] == small["strand"]
            and ov(small["i_start"], small["i_end"], big["i_start"], big["i_end"])
            >= 0.9 * (small["i_end"] - small["i_start"] + 1)
            and ov(small["j_start"], small["j_end"], big["j_start"], big["j_end"])
            >= 0.9 * (small["j_end"] - small["j_start"] + 1))


def one_arm(a, b):
    return a["strand"] == b["strand"] and (recip(a, b, 5) != recip(a, b, 3))


def key(h):
    return (h["strand"], h["i_start"], h["i_end"], h["j_start"], h["j_end"])


def main():
    fa, old_path = sys.argv[1], sys.argv[2]
    with open(old_path) as f:
        old = [dict(strand=r["Strand"], i_start=int(r["i_start"]), i_end=int(r["i_end"]),
                    j_start=int(r["j_start"]), j_end=int(r["j_end"]), score=int(r["Score"]),
                    dG=r["dG(kcal/mol)"], pp=r["percent_paired"], structure=r["structure"])
               for r in csv.DictReader(f, delimiter="\t")]
    print(f"old pipeline: {len(old)} dsRNAs")
    old_keys = {key(h): h for h in old}

    for mode, args in MODES.items():
        new, info = run_core(fa, *args)
        print(f"\n--- core {mode}: {info}")
        new_keys = {key(h): h for h in new}
        exact = old_keys.keys() & new_keys.keys()
        same_fold = sum(old_keys[k]["dG"] == new_keys[k]["dG(kcal/mol)"]
                        and old_keys[k]["structure"] == new_keys[k]["structure"] for k in exact)
        matched = [h for h in old if any(same_dsrna(h, n) for n in new)]
        print(f"  old dsRNAs with identical coordinates in core: {len(exact)}/{len(old)}; "
              f"of those, identical dG and structure: {same_fold}/{len(exact)}")
        print(f"  old dsRNAs matched (>=50% reciprocal on both arms): {len(matched)}/{len(old)} "
              f"({100 * len(matched) / len(old):.1f}%)")
        why = Counter()
        for h in old:
            if h not in matched:
                why["inside a longer core dsRNA" if any(contains(n, h) for n in new) else "other"] += 1
        for k, v in why.most_common():
            print(f"    unmatched old, {k}: {v}")
        extra = [n for n in new if not any(same_dsrna(n, h) for h in old)]
        kinds = Counter("shares one arm with an old dsRNA" if any(one_arm(n, h) for h in old)
                        else "no old dsRNA on either arm" for n in extra)
        print(f"  core dsRNAs with no old match: {len(extra)} "
              f"({', '.join(f'{k}: {v}' for k, v in kinds.most_common())})")


if __name__ == "__main__":
    main()
