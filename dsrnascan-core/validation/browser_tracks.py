"""Genome-browser tracks comparing dsRNAscan 0.5.x with dsrnascan-core on a test region.

Shifts region-relative coordinates onto the real chromosome, adds ML scores (0.5.x's
longest_helix recomputed with the fixed definition so scores are comparable), and writes
GFF3 / BEDPE / .bp files with dsRNAscan 0.5.x's own writers, for:
  old_0.5.x       the published pipeline's output
  core_default    dsrnascan-core default (reproduces 0.5.x)
  core_windowless dsrnascan-core --windowless (greed 1 + rescue)
  core_new_only   core_default dsRNAs with no 0.5.x dsRNA overlapping >=50% on both arms

Usage: python browser_tracks.py <region fasta> <old merged_results.txt> <chrom> <offset> <outdir>
  e.g. python browser_tracks.py out/chr19_300kb.fa out/old/..._merged_results.txt chr19 1000000 out/browser
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

from common import CORE, REPO

sys.path.insert(0, REPO)
from dsrnascan.dsRNAscan import (find_longest_helix, generate_bedpe_file,  # noqa: E402
                                 generate_bp_file, generate_gff3_file)
from dsrnascan.ml_scoring import score_results  # noqa: E402

COORDS = ["i_start", "i_end", "j_start", "j_end", "eff_i_start", "eff_i_end", "eff_j_start", "eff_j_end"]
STABILITY_THRESHOLD = 0.5526  # Youden threshold of the bundled stability model


def core(fa, *args):
    p = subprocess.run([CORE, fa, *args], capture_output=True, text=True, check=True)
    from io import StringIO
    return pd.read_csv(StringIO(p.stdout), sep="\t")


def matched(a, b):
    """For each row of a: does b have a dsRNA on the same strand overlapping >=50% reciprocally
    on both arms?"""
    out = np.zeros(len(a), dtype=bool)
    for s in "+-":
        bs = b[b.Strand == s]
        for n, r in a[a.Strand == s].iterrows():
            def recip(x0, x1, y0, y1):
                o = np.clip(np.minimum(x1, y1) - np.maximum(x0, y0) + 1, 0, None)
                return (o >= 0.5 * (x1 - x0 + 1)) & (o >= 0.5 * (y1 - y0 + 1))
            hit = recip(r.i_start, r.i_end, bs.i_start.values, bs.i_end.values) & \
                recip(r.j_start, r.j_end, bs.j_start.values, bs.j_end.values)
            out[a.index.get_loc(n)] = hit.any()
    return out


def write(df, name, outdir, chrom, offset):
    df = df.copy()
    for c in COORDS:
        df[c] = df[c] + offset
    df["Chromosome"] = chrom
    path = os.path.join(outdir, f"{name}_merged_results.txt")
    df.to_csv(path, sep="\t", index=False)
    generate_gff3_file(path, os.path.join(outdir, f"{name}.gff3"))
    generate_bedpe_file(path, os.path.join(outdir, f"{name}.bedpe"))
    generate_bp_file(path, os.path.join(outdir, f"{name}.forward.bp"), strand_filter="+")
    generate_bp_file(path, os.path.join(outdir, f"{name}.reverse.bp"), strand_filter="-")
    return df


def main():
    fa, old_path, chrom, offset, outdir = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
    os.makedirs(outdir, exist_ok=True)
    old = pd.read_csv(old_path, sep="\t")
    old["longest_helix"] = old["structure"].map(find_longest_helix)
    sets = {"old_0.5.x": old, "core_default": core(fa), "core_windowless": core(fa, "--windowless")}
    sets = {k: score_results(v) for k, v in sets.items()}
    sets["core_new_only"] = sets["core_default"][~matched(sets["core_default"], sets["old_0.5.x"])]

    print(f"\n{'set':16s} {'dsRNAs':>7s}  stability >= {STABILITY_THRESHOLD}   median stability")
    for name, df in sets.items():
        write(df, name, outdir, chrom, offset)
        s = df["stability_model_score"]
        print(f"{name:16s} {len(df):7d}  {(s >= STABILITY_THRESHOLD).sum():6d} ({100 * (s >= STABILITY_THRESHOLD).mean():4.1f}%)"
              f"   {s.median():.3f}")

    top = sets["core_new_only"].sort_values("stability_model_score", ascending=False).head(10)
    print("\nnew-only dsRNAs with the highest stability scores (to jump to in the browser):")
    for _, r in top.iterrows():
        print(f"  {chrom}:{r.i_start + offset:,}-{r.j_end + offset:,} ({r.Strand})  arms {r.i_end - r.i_start + 1}/"
              f"{r.j_end - r.j_start + 1} nt, dG {r['dG(kcal/mol)']}, stability {r.stability_model_score:.3f}")
    print(f"\nfiles in {os.path.abspath(outdir)}")


if __name__ == "__main__":
    main()
