"""Parse the NNDB Turner 2004 files in nndb_raw/ into turner2004.json.

Run:  python -I build_turner2004.py
Inputs are untrusted downloads: they are only read as text and parsed with
regular expressions and float().  See nndb_raw/SOURCES.md for the origin of
each file.

Pair naming everywhere in the JSON: a pair is two letters, <base on strand 1><base on
strand 2>, e.g. "GU" = G on strand 1 paired with U on strand 2.  Strand 1 is
always read 5'->3' and its partner strand 3'->5' (the drawing in NNDB).
"""
import html
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "nndb_raw"
PAR = RAW / "params"
NT = "ACGU"
PAIRS = ["AU", "CG", "GC", "GU", "UA", "UG"]
ROT = lambda p: p[::-1]  # (a,b) seen from the other strand is (b,a)


def lines(path):
    """Non-blank, non-comment lines of an NNDB text table."""
    out = []
    for ln in Path(path).read_text().splitlines():
        s = ln.strip()
        if s and not s.startswith("#"):
            out.append(s)
    return out


def num(tok):
    return None if tok == "." else float(tok) + 0.0  # + 0.0 turns -0.0 into 0.0


def html_text(path):
    t = Path(path).read_text()
    t = re.sub(r"(?s)<(script|style).*?</\1>", "", t)
    t = re.sub(r"<br\s*/?>", "\n", t)
    t = re.sub(r"</tr>", "@@ROW@@", t)
    t = re.sub(r"</td>", "|", t)
    t = html.unescape(re.sub(r"<[^>]+>", "", t))
    return t.replace("–", "-").replace("−", "-")


# ---------------------------------------------------------------- stacks
def parse_stack_html():
    """2-decimal Watson-Crick + GU stacks from the NNDB html tables."""
    stack = {p: {q: None for q in PAIRS} for p in PAIRS}
    special = {}
    for name in ("turner_2004_watson_crick_stack.htm", "turner_2004_gu_pair.htm"):
        txt = html_text(PAR / name)
        for row in txt.split("@@ROW@@")[1:]:
            cells = [c.strip() for c in row.split("|")]
            if len(cells) < 3 or not re.search(r"\d\.\d\d", cells[1]):
                continue
            seqs = re.findall(r"[ACGU]+", cells[0].replace("5'", " ").replace("3'", " "))
            if len(seqs) != 2:
                continue
            top, bot = seqs
            val = float(re.search(r"[-+]?\d+\.\d+", cells[1]).group())
            if "(-0.5)" in cells[1].replace(" ", ""):
                val = -0.5  # NNDB note (a): GG/UU set to -0.5 for prediction
            if len(top) == 4:  # 5'GGUC3'/3'CUGG5' three-stack special case
                special[top + "/" + bot] = val
                continue
            p1, p2 = top[0] + bot[0], top[1] + bot[1]
            stack[p1][p2] = val
            stack[ROT(p2)][ROT(p1)] = val  # same stack read from the other end
    return stack, special


def parse_stack_txt():
    """1-decimal table of the same stacks (used only as a cross-check)."""
    L = lines(PAR / "turner_2004_watson_crick_stack_dg.txt")
    out = {}
    i = 0
    while i < len(L):
        m = re.fullmatch(r"([ACGU])X", L[i])
        if m:
            p1 = m.group(1) + L[i + 1][0]
            rows = L[i + 3:i + 7]
            for x, r in zip(NT, rows):
                for y, v in zip(NT, r.split()[1:]):
                    if v != ".":
                        out[(p1, x + y)] = float(v)
            i += 7
        else:
            i += 1
    return out


# ----------------------------------------------------- terminal mismatch
def parse_terminal_mismatch():
    L = lines(PAR / "turner_2004_terminal_mismatch_dg.txt")
    tm, i = {}, 0
    while i < len(L):
        m = re.fullmatch(r"([ACGU])X", L[i])
        if m:
            p = m.group(1) + L[i + 1][0]
            tm[p] = {x: dict(zip(NT, map(float, r.split()[1:]))) for x, r in zip(NT, L[i + 3:i + 7])}
            i += 7
        else:
            i += 1
    return tm  # tm[pair][X][Y]; X = nt 3' of the strand-1 base, Y = nt 5' of the strand-2 base


def parse_dangles():
    L = lines(PAR / "turner_2004_dangle_dg.txt")
    d3, d5, i = {}, {}, 0
    while i < len(L):
        a, b = L[i], L[i + 1]
        vals = [num(v) for v in L[i + 3].split()]
        if len(a) == 2 and len(b) == 1 and a[1] == "X":  # 5'-P X-3' over 3'-Q-5'
            tgt, key = d3, a[0] + b
        elif len(a) == 1 and len(b) == 2 and b[1] == "X":  # 5'-P-3' over 3'-Q X-5'
            tgt, key = d5, a + b[0]
        else:
            raise ValueError("unexpected dangle block: %r %r" % (a, b))
        if any(v is not None for v in vals):
            tgt[key] = dict(zip(NT, vals))
        i += 4
    return d3, d5


# ------------------------------------------------------------ loop sizes
def parse_loop_sizes():
    bulge, internal = {}, {}
    for ln in lines(PAR / "turner_2004_loop_dg.txt"):
        f = ln.split()
        if f[0].isdigit():
            n = int(f[0])
            if f[1] != ".":
                internal[str(n)] = float(f[1])
            bulge[str(n)] = float(f[2])
    return bulge, internal


def parse_internal_misc():
    t = (PAR / "turner_2004_internal_loop_dg.txt").read_text()
    sec23 = t.split("2x3 terminal mismatch:")[1].split("Other terminal mismatch:")[0]
    other = t.split("Other terminal mismatch:")[1]
    kv = lambda txt: {k: float(v) for k, v in re.findall(r"^\s*([A-Za-z/]+)\s+(-?[\d.]+)\s*$", txt, re.M)}
    m23 = kv(sec23)  # YR/AG, RY/GA, YR/GA, GG, UU
    return {
        "asymmetry_per_nt": float(re.search(r"Asymmetry\s+([\d.]+)", t).group(1)),
        "au_gu_closure": float(re.search(r"AU/GU closure:\s+([\d.]+)", t).group(1)),
        # 2x3 loops: key = closing-pair class "/" mismatch (strand-1 nt, strand-2 nt); R = purine, Y = pyrimidine;
        # the closing pair is (strand-1 base, strand-2 base).  RY/AG is 0 in NNDB; GG and UU apply to any closing pair.
        "mismatch_2x3": {"RY/AG": 0.0, "YR/AG": m23["YR/AG"], "RY/GA": m23["RY/GA"],
                         "YR/GA": m23["YR/GA"], "GG": m23["GG"], "UU": m23["UU"]},
        # all other loops except 1x(n-1): key = (strand-1 nt)(strand-2 nt) of the first mismatch
        "mismatch_generic": kv(other),
        "mismatch_1xn": 0.0,
    }


# ----------------------------------------------------- 1x1, 2x1, 2x2 tables
def parse_int11():
    L = lines(PAR / "turner_2004_int11_dg.txt")
    out, i = {}, 0
    while i < len(L):
        if L[i] == "X":
            top, bot = L[i + 1].split(), L[i + 2].split()
            left, right = top[0] + bot[0], top[1] + bot[1]
            rows = L[i + 5:i + 9]
            out.setdefault(left, {})[right] = {x: [num(v) for v in r.split()[1:]] for x, r in zip(NT, rows)}
            i += 9
        else:
            i += 1
    return out  # int11[left][right][X] = [E for Y in ACGU]


def parse_int21():
    L = lines(PAR / "turner_2004_int21_dg.txt")
    out, i = {}, 0
    while i < len(L):
        if L[i] == "X":
            top, bot = L[i + 1].split(), L[i + 2].split()
            left, right = top[0] + bot[0], top[1] + bot[1]
            fixed = L[i + 3][1]
            rows = L[i + 5:i + 9]
            out.setdefault(left, {}).setdefault(right, {})[fixed] = {
                x: [num(v) for v in r.split()[1:]] for x, r in zip(NT, rows)}
            i += 9
        else:
            i += 1
    return out  # int21[left][right][fixed][X] = [E for Y in ACGU]


def parse_int22():
    L = lines(PAR / "turner_2004_int22_dg.txt")
    out, i = {}, 0
    while i < len(L):
        m = re.fullmatch(r"([ACGU]) X1 Y1 ([ACGU])", L[i])
        if m:
            m2 = re.fullmatch(r"([ACGU]) X2 Y2 ([ACGU])", L[i + 1])
            left, right = m.group(1) + m2.group(1), m.group(2) + m2.group(2)
            y1, y2 = L[i + 2].split(), L[i + 3].split()
            assert len(y1) == 16 and len(y2) == 16
            rows = L[i + 4:i + 20]
            tab = {}
            for r in rows:
                f = r.split()
                tab[f[0]] = [num(v) for v in f[1:]]
            out.setdefault(left, {})[right] = tab
            i += 20
        else:
            i += 1
    return out  # int22[left][right][X1X2] = [E for Y1Y2 in AA,AC,AG,AU,CA,...,UU]


# ------------------------------------------------------------------ main
def main():
    stack, special = parse_stack_html()
    chk = parse_stack_txt()
    for (p1, p2), v in chk.items():
        h = stack[p1][p2]
        assert h is not None and abs(round(h, 1) - v) < 0.051, (p1, p2, h, v)
    assert all(stack[p][q] is not None for p in PAIRS for q in PAIRS), "incomplete stack table"
    stack1 = {p: {q: chk[(p, q)] for q in PAIRS} for p in PAIRS}  # the 0.1 kcal/mol plain-text table
    for p in PAIRS:
        for q in PAIRS:
            assert stack1[p][q] == stack1[ROT(q)][ROT(p)], "stack table not symmetric"
    bulge, internal = parse_loop_sizes()
    d3, d5 = parse_dangles()
    misc = parse_internal_misc()
    out = {
        "meta": {
            "source": "NNDB RNA Turner 2004 set, https://rna.urmc.rochester.edu/NNDB/ (see nndb_raw/SOURCES.md)",
            "units": "kcal/mol at 37 C",
            "pair_naming": "two letters: base on strand 1 then base on strand 2 (strand 1 read 5'->3', strand 2 3'->5')",
            "R_kcal_per_mol_K": 1.987e-3, "T_kelvin": 310.15,
            "layout": {
                "stack[p][q], stack_1dp[p][q]": "p = pair nearer the 5' end of strand 1, q = next pair; e.g. stack.GC.CG = 5'GC3'/3'CG5'",
                "terminal_mismatch[pair][X][Y]": "helix end drawn 5'->3' on top: pair = (top base)(bottom base); X = nt 3' of the top base, Y = nt 5' of the bottom base",
                "dangle3[pair][X]": "X = nt 3' of the top base of the pair (bottom strand has no extra nt)",
                "dangle5[pair][Y]": "Y = nt 5' of the bottom base of the pair (top strand has no extra nt)",
                "bulge.init[n], interior.init[n]": "n = number of unpaired nucleotides, 1..30 (interior from 4); beyond 30 extrapolate with the log coefficient",
                "interior.mismatch_generic": "key = XY: X = loop nt next to the closing pair on its 5'->3' strand, Y = loop nt next to the partner base; applied at both closing pairs (the inner pair read from the far side, strand 2 on top); 1x(n-1) loops use mismatch_1xn",
                "int11[L][R][X][y]": "L, R = outer and inner pair; X = strand-1 loop nt; y = index (A,C,G,U) of the strand-2 loop nt",
                "int21[L][R][F][X][y]": "1x2 loop, single nt X on strand 1; strand 2 reads Y then F going away from L; y = index of Y",
                "int22[L][R][X1X2][k]": "X1,X2 = strand-1/strand-2 loop nt next to L; k = 4*idx(Y1)+idx(Y2), Y1/Y2 = strand-1/strand-2 loop nt next to R",
            },
        },
        "stack": stack,  # stack[pair_5'][pair_3'] ; pairs as above, both read along strand 1
        "stack_1dp": stack1,  # NNDB plain-text table (RNAstructure format), same stacks rounded to 0.1
        "stack_special": special,  # key "GGUC/CUGG": single value replacing the 3 stacks of that motif
        # html tables give 2 decimals; the plain-text tables are 0.1-rounded versions of the same numbers,
        # so *_1dp are those scalars rounded half away from zero (derived here, not listed in NNDB text files)
        "init": {"intermolecular": 4.09, "terminal_AU_GU_per_end": 0.45, "symmetry_self_complementary": 0.43,
                 "intermolecular_1dp": 4.1, "terminal_AU_GU_per_end_1dp": 0.5},
        "terminal_mismatch": parse_terminal_mismatch(),  # [pair][X][Y]: X = 3' neighbour of strand-1 base, Y = 5' neighbour of strand-2 base
        "dangle3": d3,  # [pair][X]: X = nt 3' of the strand-1 base of the pair
        "dangle5": d5,  # [pair][X]: X = nt 5' of the strand-2 base of the pair
        "bulge": {"init": bulge, "special_C": -0.9, "init_1_html": 3.81,
                  "log_extrapolation_coeff_RT": 1.75},
        "interior": {"init": internal, "log_extrapolation_coeff": 1.08, **misc},
        "int11": parse_int11(),  # [left pair][right pair][X][Y index in ACGU order]
        "int21": parse_int21(),  # 1 nt (X) on strand 1; strand 2 has (Y, fixed): [left][right][fixed][X][Y idx]
        "int22": parse_int22(),  # [left][right][X1X2][Y1Y2 idx in AA,AC,..,UU order]
    }
    assert sum(len(v) for v in out["int11"].values()) == 36
    assert sum(len(v) for v in out["int22"].values()) == 36
    assert sum(len(f) for v in out["int21"].values() for f in v.values()) == 144
    s = json.dumps(out, indent=1)
    s = re.sub(r"\[\s+([^\[\]]*?)\s+\]", lambda m: "[" + " ".join(m.group(1).split()) + "]", s)
    (HERE / "turner2004.json").write_text(s + "\n")
    print("wrote turner2004.json", len(s), "bytes")


if __name__ == "__main__":
    main()
