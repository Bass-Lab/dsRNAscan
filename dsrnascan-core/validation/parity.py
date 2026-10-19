"""DP parity with einverted on real sequence.

The core runs einverted's recurrence and traceback, so every hit einverted reports from a single
pass must appear among the core's unselected candidates with identical coordinates, score,
match count and gap count. Also checks that block-parallel runs give byte-identical output.

Usage: python parity.py <single-record fasta> [maxrepeat]
"""
import sys
import time
from collections import Counter

from common import COORDS, coords, read_fasta, run_core, run_einverted


def main():
    fasta = sys.argv[1]
    maxrepeat = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    seq = read_fasta(fasta)

    t = time.time()
    old = run_einverted(seq, maxrepeat)
    t_old = time.time() - t
    new, info = run_core(fasta, "--windowless", "--no-select", "--min-match-perc", 0, "--strand", "plus", "--max-span", maxrepeat,
                         "--threads", 1, "--block-size", 10**9)
    print(f"einverted single pass: {len(old)} hits in {t_old:.1f} s")
    print(f"core, 1 block, 1 thread: {len(new)} candidates | {info}")

    full = {(coords(h), h["score"], h["matches"], h["total"], h["gaps"]) for h in new}
    by_coords = {coords(h) for h in new}
    by_outer = {(h["i_start"], h["j_end"], h["score"]) for h in new}
    tally = Counter()
    misses = []
    for h in old:
        if (coords(h), h["score"], h["matches"], h["total"], h["gaps"]) in full:
            tally["identical (coords, score, matches, total, gaps)"] += 1
        elif coords(h) in by_coords:
            tally["same coords, different counts"] += 1
        elif (h["i_start"], h["j_end"], h["score"]) in by_outer:
            tally["same outer pair + score, different inner end"] += 1
            misses.append(h)
        else:
            tally["not found"] += 1
            misses.append(h)
    for k, v in tally.most_common():
        print(f"  {k}: {v}/{len(old)} ({100 * v / max(1, len(old)):.1f}%)")
    for h in misses[:10]:
        print("   e.g. einverted", {k: h[k] for k in (*COORDS, "score", "matches", "total", "gaps")})

    # Block-parallel output must be identical to the single-block run.
    ref, _ = run_core(fasta, "--windowless", "--no-select", "--min-match-perc", 0, "--strand", "plus", "--max-span", maxrepeat,
                      "--threads", 1, "--block-size", 10**9)
    for bs in (25_000, 100_000):
        par, info = run_core(fasta, "--windowless", "--no-select", "--min-match-perc", 0, "--strand", "plus", "--max-span", maxrepeat,
                             "--block-size", bs)
        same = sorted(map(str, par)) == sorted(map(str, ref))
        print(f"block size {bs:,}: {'identical' if same else 'DIFFERENT'} to single block | {info}")


if __name__ == "__main__":
    main()
