"""Shared-arm cases with known answers: one arm with a strong (exact) and a weak (10% substituted)
partner, sharing the 5' or the 3' arm, with the strong partner nearer or farther.
einverted (single pass) is expected to lose every weak pair; the core should keep all of them.

Usage: python synthetic.py
"""
import os
import random

from common import OUT, ov, run_core, run_einverted, write_fasta

COMP = str.maketrans("ACGT", "TGCA")
rng = random.Random(1)


def rc(s):
    return s.translate(COMP)[::-1]


def bg(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mut(s, frac=0.10):
    idx = set(rng.sample(range(len(s)), int(round(frac * len(s)))))
    return "".join(rng.choice([b for b in "ACGT" if b != c]) if i in idx else c for i, c in enumerate(s))


def build(parts):
    seq, where, pos = [], {}, 0
    for label, s in parts:
        if label:
            where[label] = (pos + 1, pos + len(s))
        seq.append(s)
        pos += len(s)
    return "".join(seq), where


A, A2, B = bg(200), bg(200), bg(200)
CASES = {
    "S5 shared 5' arm, strong partner nearer": (
        [(None, bg(500)), ("A", A), (None, bg(300)), ("rcA_strong", rc(A)),
         (None, bg(1500)), ("rcA_weak", mut(rc(A))), (None, bg(500))],
        [("A", "rcA_strong"), ("A", "rcA_weak")]),
    "S5eq shared 5' arm, two identical partners": (
        [(None, bg(500)), ("A", A), (None, bg(300)), ("rcA_1", rc(A)),
         (None, bg(1500)), ("rcA_2", rc(A)), (None, bg(500))],
        [("A", "rcA_1"), ("A", "rcA_2")]),
    "S3 shared 3' arm, strong partner nearer": (
        [(None, bg(500)), ("A_weak", mut(A)), (None, bg(1500)), ("A", A),
         (None, bg(300)), ("rcA", rc(A)), (None, bg(500))],
        [("A", "rcA"), ("A_weak", "rcA")]),
    "C control, two unrelated hairpins": (
        [(None, bg(500)), ("A", A), (None, bg(300)), ("rcA", rc(A)), (None, bg(2000)),
         ("B", B), (None, bg(300)), ("rcB", rc(B)), (None, bg(500))],
        [("A", "rcA"), ("B", "rcB")]),
    "S5swap shared 5' arm, weak partner nearer": (
        [(None, bg(500)), ("A", A), (None, bg(300)), ("rcA_weak", mut(rc(A))),
         (None, bg(1500)), ("rcA_strong", rc(A)), (None, bg(500))],
        [("A", "rcA_weak"), ("A", "rcA_strong")]),
    "S3swap shared 3' arm, weak partner nearer": (
        [(None, bg(500)), ("A", A), (None, bg(1500)), ("A_weak", mut(A)),
         (None, bg(300)), ("rcA", rc(A)), (None, bg(500))],
        [("A", "rcA"), ("A_weak", "rcA")]),
}


def found(hits, x, y):
    """A hit whose 5' arm covers >= 50% of region x and 3' arm >= 50% of region y."""
    for h in hits:
        if (ov(h["i_start"], h["i_end"], *x) >= 0.5 * (x[1] - x[0] + 1)
                and ov(h["j_start"], h["j_end"], *y) >= 0.5 * (y[1] - y[0] + 1)):
            return h
    return None


def main():
    n_ok = n_total = 0
    for name, (parts, expected) in CASES.items():
        seq, where = build(parts)
        fa = os.path.join(OUT, "synthetic.fa")
        write_fasta(fa, "case", seq)
        old = run_einverted(seq, 10000)
        new, _ = run_core(fa, "--windowless", "--strand", "plus", "--greed", 0)
        print(f"\n{name}  ({len(seq)} bp)  einverted: {len(old)} hits, core: {len(new)} hits")
        for a, b in expected:
            ho, hn = found(old, where[a], where[b]), found(new, where[a], where[b])
            n_total += 1
            n_ok += hn is not None
            desc = (f"{hn['i_start']}-{hn['i_end']} / {hn['j_start']}-{hn['j_end']} score {hn['score']}"
                    if hn else "")
            print(f"  {a:>6} x {b:<10}  einverted: {'found' if ho else 'MISSING':8} "
                  f"core: {'found' if hn else 'MISSING':8} {desc}")
        extra = [h for h in new if not any(found([h], where[a], where[b]) for a, b in expected)]
        for h in extra:
            print(f"  extra core hit: {h['i_start']}-{h['i_end']} / {h['j_start']}-{h['j_end']} "
                  f"score {h['score']} ({h['match_perc']:.0f}% match)")
    print(f"\ncore found {n_ok}/{n_total} expected pairs")


if __name__ == "__main__":
    main()
