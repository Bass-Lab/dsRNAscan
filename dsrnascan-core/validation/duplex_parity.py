"""The Rust duplexfold port must reproduce RNA.duplexfold exactly: energy, structure, .i and .j.
Also checks longest_helix/pairs from the traceback against the dot-bracket, and times both
(single thread each).

Usage: python duplex_parity.py <fasta> [<fasta> ...]   (pairs come from core candidates)
"""
import random
import subprocess
import sys
import time

import RNA

from common import CRATE, run_core

BIN = f"{CRATE}/target/release/duplexfold"


def longest_helix(structure):
    """Longest run of directly stacked pairs, from the dot-bracket."""
    left, right = structure.split("&")
    p1 = [i for i, ch in enumerate(left) if ch == "("]
    p2 = [j for j, ch in enumerate(right) if ch == ")"][::-1]
    best = cur = 0
    prev = None
    for a, b in zip(p1, p2):
        cur = cur + 1 if prev and a == prev[0] + 1 and b == prev[1] - 1 else 1
        best = max(best, cur)
        prev = (a, b)
    return best


def sample_pairs(fastas, rng):
    pools = {"20-60": [], "61-499": [], "500+": []}
    for fa in fastas:
        for h in run_core(fa, "--windowless", "--no-select", "--strand", "both", "--min-match-perc", 70, "--block-size", 50000)[0]:
            n = len(h["i_seq"])
            pools["20-60" if n <= 60 else "61-499" if n < 500 else "500+"].append((h["i_seq"], h["j_seq"], h["gaps"]))
    take = {"20-60": 1500, "61-499": 1500, "500+": 200}
    pairs = []
    for k, v in pools.items():
        rng.shuffle(v)
        pairs += v[:take[k]]
        print(f"  {k} nt arms: {min(len(v), take[k])} of {len(v)}")
    # constructed pairs with N in flanks and loops (rare in genomes, but must not crash or differ)
    comp = str.maketrans("ACGT", "TGCA")
    for _ in range(300):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(15, 80)))
        b = list(a.translate(comp)[::-1])
        for _ in range(rng.randint(1, 6)):
            b[rng.randrange(len(b))] = rng.choice("ACGTN")
        a = list(a)
        for _ in range(rng.randint(0, 3)):
            a[rng.randrange(len(a))] = "N"
        pairs.append(("".join(a), "".join(b), None))
    print(f"  constructed pairs with N: 300")
    return pairs


def run_bin(lines):
    t = time.time()
    p = subprocess.run([BIN], input="".join(lines), capture_output=True, text=True, check=True)
    return [line.split("\t") for line in p.stdout.rstrip("\n").split("\n")], time.time() - t


def main():
    rng = random.Random(11)
    print("sampling pairs:")
    pairs = sample_pairs(sys.argv[1:], rng)

    t = time.time()
    ref = [RNA.duplexfold(a, b) for a, b, _ in pairs]
    t_vienna = time.time() - t
    got, t_rust = run_bin(f"{a}\t{b}\n" for a, b, _ in pairs)

    n = len(pairs)
    same = {"energy": 0, "structure": 0, "i, j": 0, "all four": 0, "longest_helix": 0, "pairs": 0}
    bad = []
    for (a, b, _), r, g in zip(pairs, ref, got):
        e_ok = f"{r.energy:.2f}" == g[0]
        s_ok = r.structure == g[1]
        ij_ok = (r.i, r.j) == (int(g[2]), int(g[3]))
        same["energy"] += e_ok
        same["structure"] += s_ok
        same["i, j"] += ij_ok
        same["all four"] += e_ok and s_ok and ij_ok
        same["longest_helix"] += longest_helix(r.structure) == int(g[5])
        same["pairs"] += r.structure.count("(") == int(g[4])
        if not (e_ok and s_ok and ij_ok) and len(bad) < 5:
            bad.append((a[:40], b[:40], f"{r.energy:.2f}", g[0], r.structure[:60], g[1][:60], (r.i, r.j), g[2:4]))
    print(f"\n{n} pairs vs RNA.duplexfold (ViennaRNA {RNA.__version__}):")
    for k, v in same.items():
        print(f"  identical {k}: {v}/{n} ({100 * v / n:.2f}%)")
    for x in bad:
        print("  MISMATCH", x)
    print(f"time, single thread: ViennaRNA {t_vienna:.1f} s, Rust port (full) {t_rust:.1f} s")

    # Banded: only the anti-diagonal band around the search path, widened by a margin.
    real = [(k, p) for k, p in enumerate(pairs) if p[2] is not None]
    t_ref = sum(1 for _ in real)
    print(f"\nbanded fold on the {t_ref} real candidate pairs (vs RNA.duplexfold):")
    for margin in (0, 10, 20, 30, 60):
        got_b, t_b = run_bin(f"{a}\t{b}\t{g}\t{margin}\n" for _, (a, b, g) in real)
        ok = diff_e = 0
        worst = 0.0
        for (k, _), g in zip(real, got_b):
            r = ref[k]
            exact = f"{r.energy:.2f}" == g[0] and r.structure == g[1] and (r.i, r.j) == (int(g[2]), int(g[3]))
            ok += exact
            if not exact and g[0] != "NA":
                diff_e += f"{r.energy:.2f}" != g[0]
                worst = max(worst, float(g[0]) - r.energy)
        print(f"  margin {margin:>2}: identical {ok}/{t_ref} ({100 * ok / t_ref:.2f}%), energy differs in "
              f"{diff_e} (worst +{worst:.2f} kcal/mol), {t_b:.1f} s")


if __name__ == "__main__":
    main()
