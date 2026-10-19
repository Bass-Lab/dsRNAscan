"""ViennaRNA's Turner 2004 tables (37 C free energies and enthalpies) in src/duplex.rs's layout.

Reads the parameter file ViennaRNA writes (`RNA.params_save`, ViennaRNA 2.7.0) and maps each table
onto the layout of ../src/turner2004.rs, following duplexfold's calls into E_IntLoop/E_ExtLoop.
ViennaRNA folds the AU/GU closure into its interior-loop mismatch tables, so those are kept per
closing pair here (MMI, MM23I, MM1NI) and rescaled as whole entries, as ViennaRNA does.

The mapping is checked on the 37 C free energies: they must reproduce the validated NNDB-derived
tables (`check_37`). The same mapping then gives the enthalpies.
"""
import re

import RNA

NT = "ACGU"
PAIRS = ["AU", "CG", "GC", "GU", "UA", "UG"]  # src/duplex.rs pair index
VT = {"CG": 1, "GC": 2, "GU": 3, "UG": 4, "AU": 5, "UA": 6}  # ViennaRNA pair type
REV = lambda p: p[::-1]
B = lambda x: NT.index(x) + 1  # ViennaRNA base index


def read_par(path):
    out, cur = {}, None
    for ln in open(path):
        if ln.startswith("# "):
            cur = ln[2:].strip()
            out[cur] = []
            continue
        if cur is None:
            continue
        ln = re.sub(r"/\*.*?\*/", "", ln).strip()
        for t in ln.split():
            if t == "DEF":
                out[cur].append(-50)
            elif t == "INF":
                out[cur].append(10_000_000)
            elif re.fullmatch(r"-?\d+", t):
                out[cur].append(int(t))
            elif re.fullmatch(r"-?\d+\.\d+", t):
                out[cur].append(float(t))
    return out


def shaped(flat, dims):
    """Nested-list view of a C array stored in row-major order."""
    if len(dims) == 1:
        return list(flat[: dims[0]])
    step = len(flat) // dims[0]
    return [shaped(flat[k * step:(k + 1) * step], dims[1:]) for k in range(dims[0])]


def tables(par, suffix=""):
    """Our-layout tables from ViennaRNA's (suffix "" = 37 C, "_enthalpies" = dH)."""
    g = lambda name, dims: shaped(par[name + suffix], dims)
    stack = g("stack", [7, 7])
    mmi, mm23, mm1n = (g(n, [7, 5, 5]) for n in ("mismatch_internal", "mismatch_internal_23", "mismatch_internal_1n"))
    ext = g("mismatch_exterior", [7, 5, 5])
    d5, d3 = g("dangle5", [7, 5]), g("dangle3", [7, 5])
    i11, i21, i22 = g("int11", [7, 7, 5, 5]), g("int21", [7, 7, 5, 5, 5]), g("int22", [6, 6, 4, 4, 4, 4])
    bulge, internal = g("bulge", [31]), g("internal", [31])
    v = lambda p: VT[p] - 1  # row in the parameter file (no NP row)
    t = {
        "STACK": [[stack[v(o)][v(REV(i))] for i in PAIRS] for o in PAIRS],
        # helix ends: E_ExtLoop(type, 5' neighbour, 3' neighbour); our pair = (top, bottom) with X 3' of top
        "TMM": [[[ext[v(REV(p))][B(y)][B(x)] for y in NT] for x in NT] for p in PAIRS],
        "D3": [[d3[v(REV(p))][B(x)] for x in NT] for p in PAIRS],
        "D5": [[d5[v(REV(p))][B(y)] for y in NT] for p in PAIRS],
        "BULGE": [0] + bulge[1:31],
        "INT_INIT": [0] * 4 + internal[4:31],
        "MMI": [[[mmi[v(p)][B(x)][B(y)] for y in NT] for x in NT] for p in PAIRS],
        "MM23I": [[[mm23[v(p)][B(x)][B(y)] for y in NT] for x in NT] for p in PAIRS],
        "MM1NI": [[[mm1n[v(p)][B(x)][B(y)] for y in NT] for x in NT] for p in PAIRS],
        "INT11": [[[[i11[v(o)][v(REV(i))][B(x)][B(y)] for y in NT] for x in NT] for i in PAIRS] for o in PAIRS],
        "INT21": [[[[[i21[v(o)][v(REV(i))][B(x)][B(f)][B(y)] for y in NT] for x in NT] for f in NT]
                   for i in PAIRS] for o in PAIRS],
        "INT22": [[[[i22[v(o)][v(REV(i))][B(a) - 1][B(c) - 1][B(d) - 1][B(b) - 1]
                     for c in NT for d in NT] for a in NT for b in NT] for i in PAIRS] for o in PAIRS],
    }
    misc = par["Misc"]  # DuplexInit dG dH, TerminalAU dG dH, LXC ...
    ninio = par["NINIO"]  # m, m_dH, max
    k = 1 if suffix else 0
    t["INIT"], t["TERM_AU"], t["NINIO"] = misc[0 + k], misc[2 + k], ninio[0 + k]
    t["NINIO_MAX"] = ninio[2]  # not rescaled
    return t


def check_37(t, rust):
    """Every table of the validated 37 C set (parsed from turner2004.rs) must come out of the mapping."""
    au = lambda p: rust["INT_AU"] if p in ("AU", "GU", "UA", "UG") else 0
    want = {k: rust[k] for k in ("STACK", "TMM", "D3", "D5", "BULGE", "INT_INIT", "INT11", "INT21", "INT22",
                                 "INIT", "TERM_AU")}
    want["NINIO"], want["NINIO_MAX"] = rust["ASYM_PER_NT"], rust["ASYM_CAP"]
    want["MMI"] = [[[rust["MM_GENERIC"][x][y] + au(p) for y in range(4)] for x in range(4)] for p in PAIRS]
    want["MM23I"] = [[[rust["MM23"][k][x][y] + au(p) for y in range(4)] for x in range(4)] for k, p in enumerate(PAIRS)]
    want["MM1NI"] = [[[rust["MM_1XN"] + au(p) for y in range(4)] for x in range(4)] for p in PAIRS]
    bad = [k for k in want if t[k] != want[k]]
    assert not bad, f"ViennaRNA 37 C tables do not map onto the validated set: {bad}"


def parse_rust(path):
    src = open(path).read()
    out = {}
    for name, val in re.findall(r"pub const (\w+): [^=]+= (.+?);\n", src):
        out[name] = eval(val)
    return out


def load(par_path):
    RNA.params_load_RNA_Turner2004()
    RNA.params_save(par_path)
    return read_par(par_path)
