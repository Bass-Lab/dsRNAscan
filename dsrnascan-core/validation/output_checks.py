"""Output integrity before scaling up.

A. Every reported dsRNA is checked against the genome and against itself (coordinates, strand,
   sequences, structure, pairing, counts, energies, filters). The same checks run on dsRNAscan
   0.5.x's output as a control: a check that 0.5.x also fails is a check, not a core bug.
B. Downstream tools on the core's output: ML scoring, and 0.5.x's GFF3/BEDPE/.bp writers, diffed
   against the same tools run on 0.5.x's output.
C. dsrna-browse (browse_results.py) started headless on the core's output; its API is queried.

Usage: python output_checks.py <single-record fasta> <0.5.x merged_results.txt> <outdir> [core args, e.g. --start N --end N]
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from collections import Counter

import pandas as pd
import RNA

from common import CORE, REPO, read_fasta

sys.path.insert(0, REPO)
from dsrnascan.dsRNAscan import (find_longest_helix, generate_bedpe_file,  # noqa: E402
                                 generate_bp_file, generate_gff3_file)
from dsrnascan.ml_scoring import score_results_file  # noqa: E402

COMP = str.maketrans("ACGTN", "TGCAN")
PAIRS = {"AU", "UA", "GC", "CG", "GU", "UG"}


def rna(s):
    return s.upper().replace("T", "U")


def genome_arm(seq, start, end, strand):
    s = seq[start - 1:end]
    return rna(s if strand == "+" else s.translate(COMP)[::-1])


def check_row(r, seq):
    """Names of the checks this dsRNA fails."""
    bad = []
    i0, i1, j0, j1 = r.i_start, r.i_end, r.j_start, r.j_end
    e = (r.eff_i_start, r.eff_i_end, r.eff_j_start, r.eff_j_end)
    if not (1 <= i0 <= i1 < j0 <= j1 <= len(seq)):
        bad.append("coordinates ordered and inside the sequence")
    if not (i0 <= e[0] <= e[1] <= i1 and j0 <= e[2] <= e[3] <= j1):
        bad.append("trimmed duplex (eff_*) inside the arms")
    if r.i_seq != genome_arm(seq, e[0], e[1], r.Strand) or r.j_seq != genome_arm(seq, e[2], e[3], r.Strand):
        bad.append("i_seq/j_seq = genome at eff_* (strand-aware)")
    left, right = r.structure.split("&")
    if (len(left), len(right)) != (len(r.i_seq), len(r.j_seq)):
        bad.append("arm sequence lengths = structure halves (Forna)")
    else:
        lp = [a for a, c in enumerate(left) if c == "("]
        rp = [b for b, c in enumerate(right) if c == ")"][::-1]
        if len(lp) != len(rp) or set(left) - set("(.") or set(right) - set(")."):
            bad.append("structure balanced, intermolecular only")
        elif any(r.i_seq[a] + r.j_seq[b] not in PAIRS for a, b in zip(lp, rp)):
            bad.append("every pair canonical (AU/GC/GU)")
        bp = len(lp)
        if bp != r.base_pairs or abs(round(bp * 2 / (len(left) + len(right)) * 100, 2) - r.percent_paired) > 1e-6:
            bad.append("base_pairs and percent_paired recompute")
        if r.longest_helix != find_longest_helix(r.structure):
            bad.append("longest_helix = stacked-helix length")
        fc = RNA.fold_compound(f"{r.i_seq}&{r.j_seq}")
        if abs(fc.eval_structure(left + right) - r["dG(kcal/mol)"]) > 0.005:
            bad.append("dG = ViennaRNA energy of this structure on these sequences")
    m, t = map(int, r.RawMatch.split("/"))
    if r.Score - (3 * m - 4 * (t - m) - 12 * r.Gaps) not in (0, 3):  # +3: einverted's unreported adjacent pair
        bad.append("Score consistent with matches/mismatches/gaps")
    if abs(round(m / t * 100, 2) - r.PercMatch) > 1e-6:
        bad.append("PercMatch = matches/total")
    if not (r.PercMatch >= 70 and r.percent_paired >= 70 and r.base_pairs >= 25):
        bad.append("filters hold (match >= 70, paired >= 70, bp >= 25)")
    return bad


def check_table(df, seq, name):
    fails = Counter()
    for r in df.itertuples(index=False):
        r = pd.Series(r._asdict()).rename({"_10": "dG(kcal/mol)"})
        for b in check_row(r, seq):
            fails[b] += 1
    keys = ["Strand", "i_start", "i_end", "j_start", "j_end"]
    dup = int(df.duplicated(keys).sum())
    nested = 0
    for s, g in df.groupby("Strand"):
        a = g[keys[1:]].to_numpy()
        for x in a:
            inside = (a[:, 0] <= x[0]) & (a[:, 1] >= x[1]) & (a[:, 2] <= x[2]) & (a[:, 3] >= x[3])
            nested += int(inside.sum() > 1)
    if dup:
        fails["no duplicate dsRNAs"] += dup
    if nested:
        fails["no dsRNA nested in another (same strand)"] += nested
    print(f"\n[A] {name}: {len(df)} dsRNAs")
    names = ["coordinates ordered and inside the sequence", "trimmed duplex (eff_*) inside the arms",
             "i_seq/j_seq = genome at eff_* (strand-aware)", "arm sequence lengths = structure halves (Forna)",
             "structure balanced, intermolecular only", "every pair canonical (AU/GC/GU)",
             "base_pairs and percent_paired recompute", "longest_helix = stacked-helix length",
             "dG = ViennaRNA energy of this structure on these sequences",
             "Score consistent with matches/mismatches/gaps", "PercMatch = matches/total",
             "filters hold (match >= 70, paired >= 70, bp >= 25)", "no duplicate dsRNAs",
             "no dsRNA nested in another (same strand)"]
    for k in names:
        print(f"    {'PASS' if not fails[k] else 'FAIL':4s}  {k}" + (f"  ({fails[k]} fail)" if fails[k] else ""))
    return fails


def write_set(df, outdir, base):
    """Score with the bundled ML models and write the same files 0.5.x writes."""
    os.makedirs(outdir, exist_ok=True)
    merged = os.path.join(outdir, f"{base}_merged_results.txt")
    df.to_csv(merged, sep="\t", index=False)
    score_results_file(merged)
    generate_gff3_file(merged, os.path.join(outdir, f"{base}.gff3"))
    generate_bedpe_file(merged, os.path.join(outdir, f"{base}.bedpe"))
    generate_bp_file(merged, os.path.join(outdir, f"{base}.forward.bp"), strand_filter="+")
    generate_bp_file(merged, os.path.join(outdir, f"{base}.reverse.bp"), strand_filter="-")
    return merged


def browser_check(outdir, n_expected, port=8765):
    proc = subprocess.Popen([sys.executable, "-m", "dsrnascan.browse_results", outdir, "--no-browser", "--port", str(port)],
                            cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    get = lambda p: json.loads(urllib.request.urlopen(f"http://localhost:{port}{p}", timeout=10).read())
    try:
        for _ in range(60):
            try:
                summary = get("/api/summary")
                break
            except Exception:
                time.sleep(0.5)
        else:
            print("    FAIL  browser did not start:", proc.stdout.read()[:500] if proc.poll() is not None else "")
            return
        page = get("/api/results?page_size=200")
        ok_total = page["total"] == n_expected
        print(f"    {'PASS' if ok_total else 'FAIL'}  /api/results total = {page['total']} (expected {n_expected}); "
              f"columns: {', '.join(sorted(page['results'][0]))}")
        bad = 0
        for r in page["results"][:200]:
            d = get(f"/api/result/{r['id']}")
            left, right = d["structure"].split("&")
            if len(left) != len(d["i_seq"]) or len(right) != len(d["j_seq"]):
                bad += 1
        print(f"    {'PASS' if not bad else 'FAIL'}  /api/result/<id>: Forna sequence and structure lengths agree "
              f"for {200 - bad}/200 dsRNAs (both strands)")
        print(f"    INFO  summary keys: {', '.join(sorted(summary))}")
    finally:
        proc.terminate()
        proc.wait()


def main():
    fa, old_path, outdir, core_args = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    seq = read_fasta(fa)
    shutil.rmtree(outdir, ignore_errors=True)
    core_tsv = subprocess.run([CORE, fa, *core_args], capture_output=True, text=True, check=True).stdout
    core = pd.read_csv(pd.io.common.StringIO(core_tsv), sep="\t")
    old = pd.read_csv(old_path, sep="\t")

    check_table(core, seq, "dsrnascan-core (default)")
    check_table(old, seq, "dsRNAscan 0.5.x (control)")

    print("\n[B] ML scoring + 0.5.x writers on both outputs, files diffed")
    base = os.path.basename(old_path).replace("_merged_results.txt", "")
    old_fixed = old.copy()
    old_fixed["longest_helix"] = old_fixed["structure"].map(find_longest_helix)  # same definition as the core
    m_core = write_set(core, os.path.join(outdir, "core"), base)
    m_old = write_set(old_fixed, os.path.join(outdir, "old"), base)
    sc, so = pd.read_csv(m_core, sep="\t"), pd.read_csv(m_old, sep="\t")
    keys = ["Strand", "i_start", "i_end", "j_start", "j_end"]
    j = so.merge(sc, on=keys, suffixes=("_o", "_c"))
    for col in ("stability_model_score", "probing_model_score", "likely_edited", "likely_forms"):
        same = (j[col + "_o"].astype(str) == j[col + "_c"].astype(str)).sum()
        print(f"    {'PASS' if same == len(j) else 'FAIL'}  {col}: identical for {same}/{len(j)}")
    # 0.5.x orders dsRNAs that share an i_start by worker completion, so its sequential ID labels
    # are not reproducible even between two 0.5.x runs; compare records with the labels removed.
    unlabel = lambda line: re.sub(r"dsRNA_\S+?_\d{6}(?!\d)", "ID", line)
    for ext in (".gff3", ".bedpe", ".forward.bp", ".reverse.bp"):
        a = open(os.path.join(outdir, "old", base + ext)).read().splitlines()
        b = open(os.path.join(outdir, "core", base + ext)).read().splitlines()
        same = sorted(map(unlabel, a)) == sorted(map(unlabel, b))
        print(f"    {'PASS' if same else 'FAIL'}  {ext}: {len(b)} lines; same records as 0.5.x's: {same}; "
              f"byte-identical: {a == b}")

    print("\n[C] dsrna-browse on the core's output (headless)")
    browser_check(os.path.abspath(os.path.join(outdir, "core")), len(core))


if __name__ == "__main__":
    main()
