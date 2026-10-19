"""Validate duplex_energy.py against ViennaRNA on real dsRNAscan-core candidate duplexes.

    python validate.py [--tsv candidates.tsv] [--n 2000] [--seed 1] [--workers 8]

1. Candidates: the TSV written by dsrnascan-core (columns i_seq, j_seq); if --tsv is not given the
   core is run with  --no-select --strand both --min-match-perc 70  on validation/out/chr19_300kb.fa.
2. Sample ~n pairs stratified by arm length (short 20-60, medium 61-150 and 151-499, long 500+).
3. Per pair: RNA.duplexfold -> structure, energy, i, j; cut the strands to the structure
   (including duplexfold's one-nucleotide flanks) and compute
     (a) duplex_energy.energy  with profile "duplexfold" and with profile "nndb"
     (b) ViennaRNA's own evaluation of the same structure:
         RNA.fold_compound(t1 + "&" + t2).eval_structure(structure without "&")
   Per-loop ViennaRNA energies (fold_compound.eval_int_loop) are used to say which term is responsible
   for a difference; helix-end terms are what is left after subtracting the loops and the initiation.
4. Writes validation_results.tsv (one row per pair) and validation_report.txt next to this file.
"""
import argparse
import collections
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import RNA

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import duplex_energy as D  # noqa: E402

CORE = HERE.parent
STRATA = [("short 20-60", 20, 60, 600), ("medium 61-150", 61, 150, 550),
          ("medium 151-499", 151, 499, 550), ("long 500+", 500, 10 ** 9, 600)]
# one-at-a-time switches from the strict NNDB profile to the duplexfold value
CONVENTIONS = {
    "stacks 2-decimal -> 0.1-rounded table": {"stack": "stack_1dp"},
    "initiation 4.09 -> 4.1": {"init": D.PROFILES["duplexfold"]["init"]},
    "AU/GU end penalty 0.45 -> 0.5": {"terminal_au": D.PROFILES["duplexfold"]["terminal_au"]},
    "no special C bulge": {"special_c_bulge": False},
    "no bulge state count (-RT ln n)": {"bulge_states": False},
    "no GGUC/CUGG 3-stack motif": {"gguc_special": False},
    "asymmetry penalty capped at 3.0": {"asym_cap": 3.0},
    "interior GG/UU/2x3 mismatch values": {k: D.PROFILES["duplexfold"][k] for k in ("mm_generic", "mm_2x3", "mm_2x3_pair")},
}


PAIRS = ["AU", "UA", "CG", "GC", "GU", "UG"]


def vienna(a, b, st):
    x, y = st.split("&")
    return RNA.fold_compound(a + "&" + b).eval_structure(x + y)


def loop_case(outer, inner, l1, m2):
    """outer pair, loop (l1 on strand 1, m2 on strand 2 in 5'->3' order), inner pair; no flanks."""
    return (outer[0] + l1 + inner[0], inner[1] + m2 + outer[1], "(" + "." * len(l1) + "(&)" + "." * len(m2) + ")")


def synthetic_cases():
    """Exhaustive / random constructed duplexes, one feature each, so a mismatch names its table entry."""
    import itertools
    import random
    rng = random.Random(1)
    nts = "ACGU"
    c = {}
    # helix ends: every end pair type and every flank pattern (absent or A/C/G/U) at both ends of a 2-bp helix
    c["helix ends (22500)"] = [(f5 + p[0] + q[0] + f3, g5 + q[1] + p[1] + g3,
                                "." * len(f5) + "((" + "." * len(f3) + "&" + "." * len(g5) + "))" + "." * len(g3))
                               for p, q in itertools.product(PAIRS, PAIRS)
                               for f5, g3, f3, g5 in itertools.product(["", *nts], repeat=4)]
    for name, (n1, n2) in {"1x1 table": (1, 1), "1x2 table": (1, 2), "2x1 (rotated 1x2)": (2, 1), "2x2 table": (2, 2)}.items():
        c[name + " (every entry)"] = [loop_case(o, i, "".join(x[:n1]), "".join(x[n1:]))
                                      for o, i in itertools.product(PAIRS, PAIRS)
                                      for x in itertools.product(nts, repeat=n1 + n2)]
    for n1, n2 in [(3, 3), (2, 4), (2, 3), (3, 2), (1, 3), (1, 5)]:  # generic / 2x3 / 1xn mismatch rules, all closing pairs
        c["interior %dx%d (first mismatches x all pairs)" % (n1, n2)] = [
            loop_case(o, i, x1 + "C" * (n1 - 2) + y2 if n1 > 1 else x1, x2 + "C" * (n2 - 2) + y1 if n2 > 1 else x2)
            for o, i in itertools.product(PAIRS, PAIRS) for x1, y1, x2, y2 in itertools.product(nts, repeat=4)
            if (n1 > 1 or x1 == y2) and (n2 > 1 or y1 == x2)]
    c["bulges 1-30 nt, both strands"] = [
        loop_case(o, i, "".join(rng.choice(nts) for _ in range(n)) if side == 0 else "", "".join(rng.choice(nts) for _ in range(n)) if side else "")
        for o, i in itertools.product(PAIRS, PAIRS) for n in (1, 2, 3, 4, 6, 10, 30) for side in (0, 1) for _ in range(4)]
    c["large asymmetric / symmetric loops"] = [
        loop_case(rng.choice(PAIRS), rng.choice(PAIRS), "".join(rng.choice(nts) for _ in range(a)), "".join(rng.choice(nts) for _ in range(b)))
        for a, b in [(1, 9), (1, 28), (2, 20), (5, 12), (12, 5), (7, 7), (15, 15), (3, 26), (14, 16)] for _ in range(100)]
    n_cases = []
    for _ in range(3000):  # N in flanks and loops
        o, i = rng.choice(PAIRS), rng.choice(PAIRS)
        a, b = rng.randint(0, 4), rng.randint(0, 4)
        if a + b == 0:
            a = 1
        l1, m2 = [rng.choice("ACGUN") for _ in range(a)], [rng.choice("ACGUN") for _ in range(b)]
        s1, s2, st = loop_case(o, i, "".join(l1), "".join(m2))
        f = [rng.choice(["", "N", "A", "G"]) for _ in range(4)]
        n_cases.append((f[0] + s1 + f[1], f[2] + s2 + f[3], "." * len(f[0]) + st.split("&")[0] + "." * len(f[1]) + "&" +
                        "." * len(f[2]) + st.split("&")[1] + "." * len(f[3])))
    c["N in flanks and loops (3000)"] = n_cases
    return c


def synthetic():
    lines = ["G. constructed duplexes: profile duplexfold vs fold_compound.eval_structure (max |diff| over cases)"]
    for name, cases in synthetic_cases().items():
        diffs = [D.energy(a, b, st) - vienna(a, b, st) for a, b, st in cases]
        lines.append("   %-52s n=%-6d max|diff| %.4f  cases off by > 0.005: %d" % (
            name, len(cases), max(abs(d) for d in diffs), sum(abs(d) > 0.005 for d in diffs)))
    return lines


def work(job):
    idx, s1, s2 = job
    s1, s2 = s1.upper().replace("T", "U"), s2.upper().replace("T", "U")
    r = RNA.duplexfold(s1, s2)
    t1, t2 = D.slice_duplex(s1, s2, r.structure, r.i, r.j)
    a, b = r.structure.split("&")
    fc = RNA.fold_compound(t1 + "&" + t2)
    v_total = fc.eval_structure(a + b)
    _, _, pairs = D.parse(t1, t2, r.structure)
    n1 = len(t1)
    v_loops = [fc.eval_int_loop(k + 1, n1 + l + 1, i + 1, n1 + j + 1) / 100 for (k, l), (i, j) in zip(pairs, pairs[1:])]
    row = dict(idx=idx, len1=len(s1), len2=len(s2), n_bp=len(pairs), E_ref=r.energy, E_eval=v_total)
    for prof in ("duplexfold", "nndb"):
        tl = D.terms(t1, t2, r.structure, prof)
        row["E_mine_" + prof] = round(sum(v for _, v in tl), 2)
        mine = collections.defaultdict(float)  # per loop and for init+ends
        for lab, v in tl:
            mine[int(lab[1:lab.index("]")]) if lab.startswith("[") else -1] += v
        cat = collections.defaultdict(float)
        worst = (0.0, "")
        for t, (k, l) in enumerate(pairs[:-1]):
            i, j = pairs[t + 1]
            n_a, n_b = i - k - 1, l - j - 1
            c = "stack" if not (n_a or n_b) else "bulge" if not (n_a and n_b) else "interior"
            d = mine[t] - v_loops[t]
            cat[c] += d
            if abs(d) > abs(worst[0]) + 1e-9:
                worst = (d, "%s %dx%d %s->%s mine %.2f vienna %.2f" % (
                    c, n_a, n_b, t1[k] + t2[l], t1[i] + t2[j], mine[t], v_loops[t]))
        cat["init+ends"] = mine[-1] - (v_total - sum(v_loops))
        ends = "; ".join("%s %.2f" % (lab, v) for lab, v in tl if not lab.startswith("["))
        top = max(cat, key=lambda c: abs(cat[c]))
        row[prof + "_resp"] = top if abs(cat[top]) > 0.005 else ""
        row[prof + "_resp_detail"] = (worst[1] if top in ("stack", "bulge", "interior") else ends) if row[prof + "_resp"] else ""
        for c in ("init+ends", "stack", "bulge", "interior"):
            row["%s_d_%s" % (prof, c)] = round(cat[c], 2)
    for name, over in CONVENTIONS.items():
        row["conv: " + name] = D.energy(t1, t2, r.structure, {**D.PROFILES["nndb"], **over})
    # effect of dropping duplexfold's flanking nucleotides from the structure
    lead1, lead2 = len(a) - len(a.lstrip(".")), len(b) - len(b.lstrip("."))
    a2, b2 = a.strip("."), b.strip(".")
    row["E_mine_noflank"] = D.energy(t1[lead1:lead1 + len(a2)], t2[lead2:lead2 + len(b2)], a2 + "&" + b2)
    return row


def sample(df, n, seed):
    """~n pairs: each stratum gets its quota (scaled to n); strata with too few candidates are
    taken whole and the shortfall is made up from the other strata that have spare pairs."""
    L = df.i_seq.str.len()
    pools = [df[(L >= lo) & (L <= hi)].sample(frac=1, random_state=seed) for _, lo, hi, _ in STRATA]
    take = [min(len(p), int(round(w * n / sum(s[3] for s in STRATA)))) for p, (_, _, _, w) in zip(pools, STRATA)]
    for k in range(len(pools)):
        add = min(n - sum(take), len(pools[k]) - take[k])
        take[k] += max(add, 0)
    return pd.concat([p.iloc[:k].assign(stratum=s[0]) for p, k, s in zip(pools, take, STRATA)])


def stats(d, tols=(0.01, 0.1, 0.5)):
    d = np.round(np.asarray(d, float), 2)
    a = np.abs(d)
    return dict(n=len(d), **{"<=%g" % t: float(np.mean(a <= t + 1e-9)) for t in tols}, mad=float(a.mean()),
                max=float(a.max()), bias=float(d.mean()))


def fmt(name, s):
    return "%-44s n=%-5d  <=0.01 %6.2f%%  <=0.1 %6.2f%%  <=0.5 %6.2f%%  MAD %.4f  max %.2f  mean(mine-ref) %+.4f" % (
        name, s["n"], 100 * s["<=0.01"], 100 * s["<=0.1"], 100 * s["<=0.5"], s["mad"], s["max"], s["bias"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv")
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--no-synthetic", action="store_true", help="skip the constructed-duplex checks")
    args = ap.parse_args()
    assert RNA.cvar.temperature == 37.0 and RNA.cvar.dangles == 2 and not RNA.cvar.noGU and not RNA.cvar.noLonelyPairs, \
        "validation assumes default ViennaRNA model settings (37 C, dangles=2)"
    tsv = args.tsv
    if not tsv:
        tsv = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        with open(tsv, "w") as fh:
            subprocess.run([str(CORE / "target/release/dsrnascan-core"), str(CORE / "validation/out/chr19_300kb.fa"),
                            "--no-select", "--strand", "both", "--min-match-perc", "70"], stdout=fh, check=True)
    df = pd.read_csv(tsv, sep="\t", usecols=["i_seq", "j_seq"])
    smp = sample(df, args.n, args.seed)
    jobs = sorted(((i, r.i_seq, r.j_seq) for i, r in zip(smp.index, smp.itertuples())),
                  key=lambda j: -len(j[1]) * len(j[2]))  # longest first for load balance
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as ex:
        rows = list(ex.map(work, jobs, chunksize=4))
    res = pd.DataFrame(rows).set_index("idx").join(smp[["stratum"]]).sort_index()
    res.to_csv(HERE / "validation_results.tsv", sep="\t", float_format="%.4g")
    out = []
    P = out.append
    P("candidates: %d pairs (%s); sampled %d (seed %d); ViennaRNA %s; %.0f s" % (
        len(df), "given --tsv" if args.tsv else "fresh dsrnascan-core run", len(res), args.seed, RNA.__version__, time.time() - t0))
    P("strata: " + ", ".join("%s=%d" % (k, v) for k, v in res.stratum.value_counts().sort_index().items()))
    P("duplex size: median %d bp, max %d bp" % (res.n_bp.median(), res.n_bp.max()))
    P("")
    P("A. is ViennaRNA's eval of the same structure identical to duplexfold's energy? (checks the reference)")
    P(fmt("fold_compound.eval_structure - duplexfold", stats(res.E_eval - res.E_ref)))
    P("")
    P("B. duplex_energy.py vs RNA.duplexfold energy   (diff = mine - duplexfold, kcal/mol)")
    for prof in ("duplexfold", "nndb"):
        P(fmt("profile %-10s all" % prof, stats(res["E_mine_" + prof] - res.E_ref)))
    P("   per stratum, profile duplexfold:")
    for name, *_ in STRATA:
        g = res[res.stratum == name]
        if len(g):
            P("   " + fmt(name, stats(g.E_mine_duplexfold - g.E_ref)))
    P("   per stratum, profile nndb:")
    for name, *_ in STRATA:
        g = res[res.stratum == name]
        if len(g):
            P("   " + fmt(name, stats(g.E_mine_nndb - g.E_ref)))
    P("")
    P("C. same against ViennaRNA fold_compound.eval_structure")
    for prof in ("duplexfold", "nndb"):
        P(fmt("profile %-10s all" % prof, stats(res["E_mine_" + prof] - res.E_eval)))
    P("")
    P("D. strict-NNDB profile: effect of each convention (value under duplexfold's convention minus NNDB value)")
    P("   %-44s %8s %12s %14s" % ("convention", "pairs", "mean, all", "mean, affected"))
    for name in CONVENTIONS:
        d = np.round(res["conv: " + name] - res.E_mine_nndb, 2)
        aff = d != 0
        P("   %-44s %8d %+12.4f %+14.4f" % (name, aff.sum(), d.mean(), d[aff].mean() if aff.any() else 0.0))
    P("")
    P("E. effect of dropping duplexfold's flanking nucleotides (end mismatch/dangle terms) from the evaluation")
    P(fmt("evaluate without flanks - duplexfold", stats(res.E_mine_noflank - res.E_ref)))
    for prof in ("duplexfold", "nndb"):
        P("")
        P("F. 10 worst pairs, profile %s (diff = mine - RNA.duplexfold)" % prof)
        dd = (res["E_mine_" + prof] - res.E_ref).round(2)
        for idx in dd.abs().sort_values(ascending=False).index[:10]:
            r = res.loc[idx]
            conv = max(CONVENTIONS, key=lambda c: abs(r["conv: " + c] - r.E_mine_nndb))
            P("   row %-6d %-15s arms %4d/%-4d bp %-4d duplexfold %9.2f mine %9.2f diff %+7.2f  term: %s %s%s" % (
                idx, r.stratum, r.len1, r.len2, r.n_bp, r.E_ref, r["E_mine_" + prof], dd[idx],
                r[prof + "_resp"] or "none (every term agrees to 0.005)", r[prof + "_resp_detail"] or "",
                ("  [largest convention: %s]" % conv) if prof == "nndb" else ""))
    if not args.no_synthetic:
        P("")
        out.extend(synthetic())
    text = "\n".join(out)
    print(text)
    (HERE / "validation_report.txt").write_text(text + "\n")


if __name__ == "__main__":
    main()
