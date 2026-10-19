"""Sweep the greed dial against dsRNAscan 0.5.x on one region: speed, volume, agreement, and how
many dsRNAs the bundled stability model scores as likely edited.

Usage: python greed_sweep.py <fa> <old *_merged_results.txt> [greed ...]
"""
import csv
import resource
import subprocess
import sys
import time
from io import StringIO

import pandas as pd

from common import CORE, REPO
from vs_pipeline import key, same_dsrna

sys.path.insert(0, REPO)
from dsrnascan.dsRNAscan import find_longest_helix  # noqa: E402
from dsrnascan.ml_scoring import score_results  # noqa: E402

THR = 0.5526


def rows(df):
    return [dict(strand=r.Strand, i_start=r.i_start, i_end=r.i_end, j_start=r.j_start, j_end=r.j_end)
            for r in df.itertuples()]


def main():
    fa, old_path = sys.argv[1], sys.argv[2]
    greeds = [float(g) for g in sys.argv[3:]] or [1.0, 0.9, 0.8, 0.7, 0.5, 0.0]
    old_df = pd.read_csv(old_path, sep="\t")
    old_df["longest_helix"] = old_df["structure"].map(find_longest_helix)
    old_pos = (score_results(old_df)["stability_model_score"] >= THR).sum()
    old = rows(old_df)
    old_keys = {key(h) for h in old}
    print(f"\n0.5.x: {len(old)} dsRNAs, {old_pos} above stability {THR}")
    print(f"{'mode':>28s} {'dsRNAs':>7s} {'wall s':>7s} {'CPU s':>7s} {'0.5.x recalled':>15s} {'identical':>9s} "
          f"{'also in 0.5.x':>14s} {'stability+':>10s}")
    for g in greeds:
        for extra in ([], ["--banded"], ["--no-rescue", "--banded"]):
            r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
            t = time.time()
            out = subprocess.run([CORE, fa, "--windowless", "--greed", str(g), *extra], capture_output=True, text=True, check=True).stdout
            wall = time.time() - t
            r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
            cpu = (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime)
            df = pd.read_csv(StringIO(out), sep="\t")
            new = rows(df)
            rec = sum(any(same_dsrna(h, n) for n in new) for h in old)
            prec = sum(any(same_dsrna(n, h) for h in old) for n in new)
            ident = len(old_keys & {key(h) for h in new})
            pos = (score_results(df)["stability_model_score"] >= THR).sum() if len(df) else 0
            mode = f"greed {g}" + "".join(" " + e.lstrip("-") for e in extra)
            print(f"{mode:>28s} {len(new):7d} {wall:7.1f} {cpu:7.1f} {rec:6d} ({100 * rec / len(old):4.1f}%) {ident:9d} "
                  f"{prec:6d} ({100 * prec / max(1, len(new)):4.1f}%) {pos:10d}", flush=True)


if __name__ == "__main__":
    main()
