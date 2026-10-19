"""Lock-in test for the einverted replay: on a sequence exactly one window long, the core
(-s, --no-select: einverted's raw window hits) must report exactly what einverted reports,
hit for hit (coordinates, score, matches, total, gaps).

Random windows from the given FASTAs, both strands (minus = complement, not reversed, which is
what dsRNAscan 0.5.x hands einverted).

Usage: python single_window.py <window size> <n windows per fasta> <fasta> [<fasta> ...]
"""
import os
import random
import sys

from common import COORDS, OUT, read_fasta, run_core, run_einverted, write_fasta

COMP = str.maketrans("ACGTN", "TGCAN")
KEY = (*COORDS, "score", "matches", "total", "gaps")


def main():
    w, n = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(5)
    fa = os.path.join(OUT, "window.fa")
    tot_e = tot_c = same = windows_ok = windows = 0
    shown = 0
    for path in sys.argv[3:]:
        seq = read_fasta(path)
        for _ in range(n):
            s = rng.randrange(0, len(seq) - w + 1)
            for strand in "+-":
                win = seq[s:s + w] if strand == "+" else seq[s:s + w].translate(COMP)
                write_fasta(fa, "w", win)
                e = {tuple(h[k] for k in KEY) for h in run_einverted(win, w)}
                c = {tuple(h[k] for k in KEY) for h in
                     run_core(fa, "--no-select", "-s", 500, "--strand", "plus", "--min-match-perc", 0, "-w", w)[0]}
                windows += 1
                tot_e += len(e)
                tot_c += len(c)
                same += len(e & c)
                windows_ok += e == c
                if e != c and shown < 5:
                    shown += 1
                    print(f"  {os.path.basename(path)} {s}-{s + w} ({strand}): einverted {len(e)}, core {len(c)}")
                    print(f"     only einverted: {sorted(e - c)[:3]}")
                    print(f"     only core:      {sorted(c - e)[:3]}")
    print(f"\n{windows} windows of {w} nt: identical hit sets in {windows_ok}/{windows}")
    print(f"hits: einverted {tot_e}, core {tot_c}, identical {same} "
          f"({100 * same / max(1, tot_e):.2f}% of einverted's)")


if __name__ == "__main__":
    main()
