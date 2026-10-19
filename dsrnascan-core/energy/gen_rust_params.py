"""Generate ../src/turner2004.rs: Turner 2004 free energies at 37 C and enthalpies.

The 37 C tables are filled by duplex_energy.py's lookup rules on the NNDB tables (profile
"duplexfold", validated against RNA.duplexfold), so they cannot drift from the reference evaluator.
The enthalpies come from ViennaRNA's parameter set (vienna_tables.py), mapped onto the same layout;
the mapping is checked by requiring ViennaRNA's 37 C values to reproduce the NNDB-derived tables.
Energies in dcal/mol (integer kcal/mol x 100), as in ViennaRNA. Run: python gen_rust_params.py
(needs the ViennaRNA Python package).
"""
import tempfile
from pathlib import Path

import duplex_energy as de
import vienna_tables as vt

T, P = de.T, de.PROFILES["duplexfold"]
NT = "ACGU"
PAIRS = vt.PAIRS


def dc(x):
    v = round(x * 100)
    assert abs(v - x * 100) < 1e-6, x
    return v


def nndb_37():
    """The validated 37 C set, in the layout src/duplex.rs used before enthalpies were added."""
    stack, tm = T[P["stack"]], T["terminal_mismatch"]
    g = {**T["interior"]["mismatch_generic"], **P["mm_generic"]}
    return {
        "INIT": dc(P["init"]), "TERM_AU": dc(P["terminal_au"]), "INT_AU": dc(T["interior"]["au_gu_closure"]),
        "ASYM_PER_NT": dc(T["interior"]["asymmetry_per_nt"]), "ASYM_CAP": dc(P["asym_cap"]),
        "MM_1XN": dc(T["interior"]["mismatch_1xn"]),
        "STACK": [[dc(stack[p][q]) for q in PAIRS] for p in PAIRS],
        "TMM": [[[dc(tm[p][x][y]) for y in NT] for x in NT] for p in PAIRS],
        "D3": [[dc(T["dangle3"][p][x]) for x in NT] for p in PAIRS],
        "D5": [[dc(T["dangle5"][p][y]) for y in NT] for p in PAIRS],
        "BULGE": [0] + [dc(T["bulge"]["init"][str(n)]) for n in range(1, 31)],
        "INT_INIT": [0] * 4 + [dc(T["interior"]["init"][str(n)]) for n in range(4, 31)],
        "MM_GENERIC": [[dc(g.get(x + y, 0.0)) for y in NT] for x in NT],
        "MM23": [[[dc(de.mm23(p, x, y, P)) for y in NT] for x in NT] for p in PAIRS],
        "INT11": [[[[dc(T["int11"][o][i][x][yi]) for yi in range(4)] for x in NT] for i in PAIRS] for o in PAIRS],
        "INT21": [[[[[dc(T["int21"][o][i][f][x][yi]) for yi in range(4)] for x in NT] for f in NT] for i in PAIRS]
                  for o in PAIRS],
        "INT22": [[[[dc(v) for v in T["int22"][o][i][x1 + y1]] for x1 in NT for y1 in NT] for i in PAIRS]
                  for o in PAIRS],
    }


def with_n(table, axes):
    """Add N (index 4) on the base axes: the maximum over A/C/G/U, as ViennaRNA fills its N entries
    (separately for free energies and enthalpies, so N rescales with its own pair of values)."""
    import itertools

    def get(t, idx):
        for i in idx:
            t = t[i]
        return t

    def dims(t):
        return [len(t)] + dims(t[0]) if isinstance(t, list) else []

    shape = [5 if k in axes else n for k, n in enumerate(dims(table))]

    def build(prefix):
        k = len(prefix)
        if k == len(shape):
            subs = [[i] if (j not in axes or i < 4) else range(4) for j, i in enumerate(prefix)]
            return max(get(table, c) for c in itertools.product(*subs))
        return [build(prefix + [i]) for i in range(shape[k])]

    return build([])


def int22_5(t):
    """INT22[o][i][4a+b][4c+d] -> [o][i][5a+b][5c+d] with N entries."""
    t4 = [[[[[[t[o][i][4 * a + b][4 * c + d] for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
           for i in range(6)] for o in range(6)]
    t5 = with_n(t4, {2, 3, 4, 5})
    return [[[[t5[o][i][a][b][c][d] for c in range(5) for d in range(5)] for a in range(5) for b in range(5)]
             for i in range(6)] for o in range(6)]


def add_n(t):
    axes = {"TMM": {1, 2}, "D3": {1}, "D5": {1}, "MMI": {1, 2}, "MM23I": {1, 2}, "MM1NI": {1, 2},
            "INT11": {2, 3}, "INT21": {2, 3, 4}}
    out = dict(t)
    for k, a in axes.items():
        out[k] = with_n(t[k], a)
    out["INT22"] = int22_5(t["INT22"])
    return out


def arr(x):
    if not isinstance(x, list):
        return str(x)
    return "[" + ", ".join(arr(v) for v in x) + "]"


def rust_type(x):
    return "i32" if not isinstance(x, list) else f"[{rust_type(x[0])}; {len(x)}]"


DOCS = {
    "INIT": "intermolecular initiation",
    "TERM_AU": "AU/GU helix-end penalty (helix ends and bulges > 1 nt)",
    "NINIO": "interior-loop asymmetry per nt",
    "STACK": "[outer pair][next pair], outer = nearer the 5' end of strand 1",
    "TMM": "terminal mismatch [pair top|bottom][X 3' of top][Y 5' of bottom]",
    "D3": "3' dangle [pair][X 3' of top]",
    "D5": "5' dangle [pair][Y 5' of bottom]",
    "BULGE": "bulge initiation by size (index 0 unused)",
    "INT_INIT": "generic interior-loop initiation by total size (indices < 4 unused)",
    "MMI": "generic interior first mismatch incl. AU/GU closure [closing pair][X on its 5'->3' strand][Y]",
    "MM23I": "2x3 interior first mismatch incl. AU/GU closure [closing pair][X][Y]",
    "MM1NI": "1xn interior first mismatch incl. AU/GU closure [closing pair][X][Y]",
    "INT11": "1x1 loops [outer][inner][strand-1 nt][strand-2 nt]",
    "INT21": "1x2 loops [outer][inner][F][X][Y]: X single nt on strand 1, strand 2 reads Y then F",
    "INT22": "2x2 loops [outer][inner][5*X1+X2 next to outer][5*Y1+Y2 next to inner]",
}

with tempfile.TemporaryDirectory() as d:
    par = vt.load(f"{d}/rna_turner2004.par")
g37, dh = vt.tables(par), vt.tables(par, "_enthalpies")
vt.check_37(g37, nndb_37())
g37, dh = add_n(g37), add_n(dh)

out = [
    "// @generated by energy/gen_rust_params.py. Do not edit. Turner 2004 (NNDB; Mathews et al. 2004 PNAS):\n",
    "// free energies at 37 C (NNDB, with RNA.duplexfold's conventions, see energy/README.md) and\n",
    "// enthalpies (ViennaRNA 2.7.0's Turner 2004 set; Lorenz et al. 2011), in dcal/mol.\n",
    "// Pair index: AU CG GC GU UA UG; base index: A C G U N (N = the worst of A/C/G/U, as in ViennaRNA).\n\n",
    f"/// interior-loop asymmetry cap (not rescaled with temperature)\npub const NINIO_MAX: i32 = {g37['NINIO_MAX']};\n",
]
for name, doc in DOCS.items():
    ty = rust_type(g37[name])
    out.append(f"/// {doc}\npub const {name}: {ty} = {arr(g37[name])};\n")
    out.append(f"pub const {name}_DH: {ty} = {arr(dh[name])};\n")

dest = Path(__file__).resolve().parent.parent / "src" / "turner2004.rs"
dest.write_text("".join(out))
print(f"wrote {dest}")
