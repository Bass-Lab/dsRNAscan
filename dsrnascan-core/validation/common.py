"""Shared helpers: run einverted and dsrnascan-core, parse both into the same hit dicts.

Hit dict keys: i_start, i_end (5' arm), j_start, j_end (3' arm), score, matches, total, gaps.
Coordinates are 1-based inclusive, as einverted prints them.
"""
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
CRATE = os.path.dirname(HERE)
REPO = os.path.dirname(CRATE)
CORE = os.path.join(CRATE, "target", "release", "dsrnascan-core")
EINVERTED = os.path.join(REPO, "dsrnascan", "tools", "platform_binaries", "einverted_macos_arm64")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

# dsRNAscan 0.5.x defaults
GAP, THRESHOLD, MATCH, MISMATCH = 12, 75, 3, -4

HDR = re.compile(r"^(\S+): Score (\d+): (\d+)/(\d+) \(\s*\d+%\) matches, (\d+) gaps")
COORDS = ("i_start", "i_end", "j_start", "j_end")


def run_einverted(seq, maxrepeat, name="s"):
    cmd = [EINVERTED, "-sequence", "stdin", "-gap", str(GAP), "-threshold", str(THRESHOLD),
           "-match", str(MATCH), "-mismatch", str(MISMATCH), "-maxrepeat", str(maxrepeat),
           "-outfile", "stdout", "-outseq", "/dev/null"]
    p = subprocess.run(cmd, input=f">{name}\n{seq}\n", capture_output=True, text=True, check=True)
    lines, hits, n = p.stdout.split("\n"), [], 0
    while n < len(lines):
        m = HDR.match(lines[n])
        if m and n + 3 < len(lines):
            a, b = lines[n + 1].split(), lines[n + 3].split()
            hits.append(dict(i_start=int(a[0]), i_end=int(a[-1]), j_start=int(b[-1]), j_end=int(b[0]),
                             score=int(m.group(2)), matches=int(m.group(3)), total=int(m.group(4)),
                             gaps=int(m.group(5))))
            n += 4
        else:
            n += 1
    return hits


def run_core(fasta, *args):
    """Returns (hits, stderr summary line). Reads both output formats: raw candidates
    (--no-select) and dsRNAs (merged_results columns, mapped to the same keys)."""
    p = subprocess.run([CORE, fasta, *map(str, args)], capture_output=True, text=True, check=True)
    lines = p.stdout.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    hits = []
    for line in lines[1:]:
        r = dict(zip(head, line.split("\t")))
        if "Score" in r:  # merged_results format
            m, t = r["RawMatch"].split("/")
            r.update(strand=r["Strand"], score=r["Score"], matches=m, total=t, gaps=r["Gaps"],
                     match_perc=r["PercMatch"])
        for k in ("i_start", "i_end", "j_start", "j_end", "score", "matches", "total", "gaps"):
            r[k] = int(r[k])
        r["match_perc"] = float(r["match_perc"])
        hits.append(r)
    return hits, p.stderr.strip()


def write_fasta(path, name, seq):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for k in range(0, len(seq), 60):
            f.write(seq[k:k + 60] + "\n")


def read_fasta(path):
    with open(path) as f:
        return "".join(l.strip() for l in f if not l.startswith(">")).upper()


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def coords(h):
    return tuple(h[k] for k in COORDS)
