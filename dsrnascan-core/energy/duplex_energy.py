"""Reference evaluator: Turner 2004 free energy of a two-strand RNA duplex structure.

    energy(s1, s2, structure) -> kcal/mol at 37 C

s1, s2      the two strands, 5'->3' (DNA letters accepted, T is read as U)
structure   "dotbracket1&dotbracket2" as returned by RNA.duplexfold: '(' on strand 1,
            ')' on strand 2, '.' unpaired; the k-th '(' pairs with the k-th ')' counted from
            the end of strand 2.  No intramolecular pairs, no multiloops.  len(part 1) must
            equal len(s1) and len(part 2) len(s2): pass the strands cut to the structure
            (slice_duplex() does this), INCLUDING the one unpaired flanking nucleotide
            duplexfold keeps on each side -- those flanks are the dangling ends / terminal
            mismatches of the two helix ends.  A flank that does not exist (strand end) is
            simply absent.

Model (NNDB Turner 2004, tables in turner2004.json):
    E = intermolecular initiation
      + each helix end: terminal mismatch if both flanking nucleotides exist, else the 3' or
        5' dangle that exists, plus the AU/GU end penalty if the end pair is AU/GU
      + nearest-neighbour stacks, bulge loops and interior loops between consecutive pairs
Pair names are "<strand-1 base><strand-2 base>", see build_turner2004.py for table layouts.

Profiles (see README):
    "duplexfold"  default.  Reproduces RNA.duplexfold / fold_compound.eval_structure.
    "nndb"        the NNDB pages taken literally (2-decimal stacks, special C bulge, bulge state
                  count, GGUC motif).  Differs from ViennaRNA by the terms listed in README.
Non-ACGU letters (N) in unpaired positions are scored as the least favourable of A/C/G/U for
the term they sit in, which is what ViennaRNA does.
terms() returns the itemised breakdown that energy() sums.
"""
import itertools
import json
import math
from pathlib import Path

T = json.loads((Path(__file__).resolve().parent / "turner2004.json").read_text())
RT = T["meta"]["R_kcal_per_mol_K"] * T["meta"]["T_kelvin"]
NT = "ACGU"
AUGU = {"AU", "UA", "GU", "UG"}
CANON = AUGU | {"CG", "GC"}
PURINE = "AG"

PROFILES = {
    "nndb": dict(stack="stack", init=T["init"]["intermolecular"], terminal_au=T["init"]["terminal_AU_GU_per_end"],
                 special_c_bulge=True, bulge_states=True, gguc_special=True, asym_cap=None,
                 mm_generic={}, mm_2x3={}, mm_2x3_pair={}),
    # 0.1-rounded NNDB plain-text tables plus the conventions measured in RNA.duplexfold.  Every entry
    # that is not an NNDB value was found by black-box comparison with RNA.duplexfold /
    # fold_compound.eval_structure output (validate.py, README), not copied from a ViennaRNA file.
    "duplexfold": dict(stack="stack_1dp", init=T["init"]["intermolecular_1dp"],
                       terminal_au=T["init"]["terminal_AU_GU_per_end_1dp"],
                       special_c_bulge=False, bulge_states=False, gguc_special=False,
                       asym_cap=3.0,                                # asymmetry penalty capped at 3.0
                       mm_generic={"GG": -1.0, "UU": -0.6},         # NNDB text: -1.2, -0.7
                       mm_2x3={"GG": -0.7, "UU": -0.3},             # NNDB text: -0.8, -0.4
                       mm_2x3_pair={("AU", "GA"): -1.1, ("GU", "GA"): -1.1}),  # NNDB RY/GA: -1.2
}


def clean(s):
    s = s.upper().replace("T", "U")
    return "".join(c if c in NT else "N" for c in s)


def parse(s1, s2, structure):
    s1, s2 = clean(s1), clean(s2)
    try:
        d1, d2 = structure.split("&")
    except ValueError:
        raise ValueError("structure must be 'dotbracket1&dotbracket2'")
    if len(d1) != len(s1) or len(d2) != len(s2):
        raise ValueError("structure parts (%d,%d) do not match sequence lengths (%d,%d)"
                         % (len(d1), len(d2), len(s1), len(s2)))
    if set(d1) - set("(.") or set(d2) - set(")."):
        raise ValueError("only '(' on strand 1 and ')' on strand 2 are allowed (no intramolecular pairs)")
    p1 = [i for i, c in enumerate(d1) if c == "("]
    p2 = [j for j, c in enumerate(d2) if c == ")"][::-1]
    if len(p1) != len(p2) or not p1:
        raise ValueError("unbalanced or empty duplex")
    pairs = list(zip(p1, p2))
    for i, j in pairs:
        if s1[i] + s2[j] not in CANON:
            raise ValueError("non-canonical pair %s%d-%s%d" % (s1[i], i + 1, s2[j], j + 1))
    return s1, s2, pairs


def slice_duplex(s1, s2, structure, i, j):
    """Cut s1, s2 to the part RNA.duplexfold's structure describes.  i, j are duplexfold's result.i
    and result.j: 1-based index of the last s1 nucleotide covered and of the first s2 nucleotide."""
    d1, d2 = structure.split("&")
    return s1[i - len(d1):i], s2[j - 1:j - 1 + len(d2)]


def nb(s, k):
    """Neighbour at index k, or None if off the strand."""
    return s[k] if 0 <= k < len(s) else None


def worst_case(fn, s1, s2, idx1, idx2):
    """fn(s1, s2) -> list of (label, kcal).  Any N at the given unpaired indices is replaced by each of
    A/C/G/U (jointly) and the highest-energy outcome is kept."""
    spots = [(0, i) for i in idx1 if s1[i] == "N"] + [(1, j) for j in idx2 if s2[j] == "N"]
    if not spots:
        return fn(s1, s2)
    best = None
    for combo in itertools.product(NT, repeat=len(spots)):
        t = [list(s1), list(s2)]
        for (w, k), c in zip(spots, combo):
            t[w][k] = c
        r = fn("".join(t[0]), "".join(t[1]))
        if best is None or sum(v for _, v in r) > sum(v for _, v in best):
            best = r
    return best


def end_terms(top, bot, X, Y, P, label):
    """Helix end drawn as in NNDB: base `top` (on the 5'->3' strand) paired with `bot`; X is the
    nucleotide 3' of `top`, Y the nucleotide 5' of `bot` (None if absent)."""
    pair = top + bot
    out = []
    if X and Y:
        out.append((label + " terminal mismatch " + pair + " " + X + Y, T["terminal_mismatch"][pair][X][Y]))
    elif X:
        out.append((label + " 3' dangle " + pair + " " + X, T["dangle3"][pair][X]))
    elif Y:
        out.append((label + " 5' dangle " + pair + " " + Y, T["dangle5"][pair][Y]))
    if pair in AUGU:
        out.append((label + " AU/GU end " + pair, P["terminal_au"]))
    return out


def init_n(table, n, coeff):
    """Loop initiation for n unpaired nt: table to 30, then logarithmic extrapolation (truncated to
    0.01 kcal/mol, as duplexfold's integer arithmetic does; NNDB does not specify the rounding)."""
    if str(n) in table:
        return table[str(n)]
    return table["30"] + math.trunc(100 * coeff * math.log(n / 30.0)) / 100


def bulge_terms(s1, s2, k, l, i, j, n1, n2, P, extent):
    n = n1 + n2
    outer, inner = s1[k] + s2[l], s1[i] + s2[j]
    out = [("bulge init %d" % n, init_n(T["bulge"]["init"], n, T["bulge"]["log_extrapolation_coeff_RT"] * RT))]
    if n == 1:
        # the helix stacks across a single-nucleotide bulge: no AU/GU end penalty
        out.append(("bulge-1 stack %s/%s" % (outer, inner), T[P["stack"]][outer][inner]))
        seq, b, lo, hi = (s1, k + 1, k, i) if n1 else (s2, j + 1, j, l)  # strand, bulged index, paired neighbours
        x = seq[b]
        if P["special_c_bulge"] and x == "C" and "C" in (seq[lo], seq[hi]):
            out.append(("special C bulge", T["bulge"]["special_C"]))
        if P["bulge_states"]:  # identical nucleotides the bulge could slide over, within the helix
            a = z = b
            while a - 1 >= extent[0] and seq[a - 1] == x:
                a -= 1
            while z + 1 <= extent[1] and seq[z + 1] == x:
                z += 1
            if z > a:
                out.append(("bulge-1 states x%d" % (z - a + 1), -RT * math.log(z - a + 1)))
    else:
        for name, pr in (("outer", outer), ("inner", inner)):
            if pr in AUGU:
                out.append(("bulge AU/GU closure " + name + " " + pr, P["terminal_au"]))
    return out


def mm23(close, x, y, P):
    """2x3 first-mismatch term; close = (strand-1 base, strand-2 base) of the closing pair."""
    m = {**T["interior"]["mismatch_2x3"], **P["mm_2x3"]}
    if (close, x + y) in P["mm_2x3_pair"]:
        return P["mm_2x3_pair"][(close, x + y)]
    if x + y in ("GG", "UU"):
        return m[x + y]
    if x + y in ("AG", "GA"):
        cls = ("R" if close[0] in PURINE else "Y") + ("R" if close[1] in PURINE else "Y")
        return m[cls + "/" + x + y]
    return 0.0


def interior_terms(s1, s2, k, l, i, j, n1, n2, P):
    outer, inner = s1[k] + s2[l], s1[i] + s2[j]
    rinner = inner[::-1]  # inner pair seen from the far side of the loop (strand 2 as 'top')
    x1, y1 = s1[k + 1], s2[l - 1]  # first mismatch next to the outer pair (loop nt on strand 1 / strand 2)
    x2, y2 = s2[j + 1], s1[i - 1]  # next to the inner pair, read from the other side (strand 2 / strand 1)
    I = T["interior"]
    name = "interior %dx%d" % (n1, n2)
    if (n1, n2) == (1, 1):
        return [(name, T["int11"][outer][inner][x1][NT.index(y1)])]
    if (n1, n2) == (1, 2):  # 1 nt on strand 1; strand 2 reads (y1, then `fixed`) away from the outer pair
        return [(name, T["int21"][outer][inner][s2[l - 2]][x1][NT.index(y1)])]
    if (n1, n2) == (2, 1):  # rotate 180 degrees so the single nucleotide is on the 'top' strand again
        return [(name, T["int21"][rinner][outer[::-1]][s1[i - 2]][x2][NT.index(y2)])]
    if (n1, n2) == (2, 2):
        return [(name, T["int22"][outer][inner][x1 + y1][NT.index(y2) * 4 + NT.index(x2)])]
    n = n1 + n2
    out = [(name + " init", init_n(I["init"], n, I["log_extrapolation_coeff"]))]
    asym = I["asymmetry_per_nt"] * abs(n1 - n2)
    if P["asym_cap"] is not None:
        asym = min(asym, P["asym_cap"])
    out.append((name + " asymmetry", asym))
    if min(n1, n2) == 1:  # 1 x (n-1): no mismatch bonus
        m1 = m2 = I["mismatch_1xn"]
    elif {n1, n2} == {2, 3}:
        m1, m2 = mm23(outer, x1, y1, P), mm23(rinner, x2, y2, P)
    else:
        g = {**I["mismatch_generic"], **P["mm_generic"]}
        m1, m2 = g.get(x1 + y1, 0.0), g.get(x2 + y2, 0.0)
    out.append((name + " mismatches", m1 + m2))
    for pr in (outer, inner):
        if pr in AUGU:
            out.append((name + " AU/GU closure " + pr, I["au_gu_closure"]))
    return out


def terms(s1, s2, structure, profile="duplexfold", symmetry=False):
    """Itemised energy: list of (label, kcal/mol)."""
    P = PROFILES[profile] if isinstance(profile, str) else profile
    s1, s2, pairs = parse(s1, s2, structure)
    out = [("intermolecular initiation", P["init"])]
    (i0, j0), (iL, jL) = pairs[0], pairs[-1]
    # 5' end of strand 1 = 3' end of strand 2: drawn with strand 2 on top
    out += worst_case(lambda a, b: end_terms(b[j0], a[i0], nb(b, j0 + 1), nb(a, i0 - 1), P, "5'(s1) end"),
                      s1, s2, [i0 - 1] if i0 else [], [j0 + 1] if j0 + 1 < len(s2) else [])
    out += worst_case(lambda a, b: end_terms(a[iL], b[jL], nb(a, iL + 1), nb(b, jL - 1), P, "3'(s1) end"),
                      s1, s2, [iL + 1] if iL + 1 < len(s1) else [], [jL - 1] if jL else [])
    stacks = {}  # index of the outer pair -> position in `out`, for the GGUC motif
    for t in range(len(pairs) - 1):  # loop t lies between pair t and pair t+1; its terms are labelled "[t] ..."
        (k, l), (i, j) = pairs[t], pairs[t + 1]
        n1, n2 = i - k - 1, l - j - 1
        if n1 == 0 and n2 == 0:
            p, q = s1[k] + s2[l], s1[i] + s2[j]
            lt = [("stack %s/%s" % (p, q), T[P["stack"]][p][q])]
            stacks[t] = len(out)
        elif n1 == 0 or n2 == 0:
            lt = worst_case(lambda a, b: bulge_terms(a, b, k, l, i, j, n1, n2, P, (i0, iL) if n1 else (jL, j0)),
                            s1, s2, range(k + 1, i), range(j + 1, l))
        else:
            lt = worst_case(lambda a, b: interior_terms(a, b, k, l, i, j, n1, n2, P),
                            s1, s2, range(k + 1, i), range(j + 1, l))
        out += [("[%d] %s" % (t, lab), v) for lab, v in lt]
    if P["gguc_special"]:  # 5'GGUC3'/3'CUGG5' = pairs GC,GU,UG,CG: one value replaces its three stacks
        for t in range(len(pairs) - 3):
            if all(t + m in stacks for m in range(3)) and \
                    [s1[pairs[t + m][0]] + s2[pairs[t + m][1]] for m in range(4)] == ["GC", "GU", "UG", "CG"]:
                for m in range(3):
                    out[stacks[t + m]] = (out[stacks[t + m]][0] + " [GGUC motif]", 0.0)
                out[stacks[t]] = ("GGUC/CUGG motif (3 stacks)", T["stack_special"]["GGUC/CUGG"])
    if symmetry and s1 == s2 and len(pairs) == len(s1) and s1 == s1[::-1].translate(str.maketrans("ACGU", "UGCA")):
        out.append(("symmetry correction", T["init"]["symmetry_self_complementary"]))
    return out


def energy(s1, s2, structure, profile="duplexfold", symmetry=False):
    """Duplex free energy in kcal/mol.  symmetry=True adds the NNDB self-complementary-duplex term
    (duplexfold never does)."""
    return round(sum(v for _, v in terms(s1, s2, structure, profile, symmetry)), 2)


if __name__ == "__main__":
    # worked examples from the NNDB Watson-Crick and GU pages (2-decimal values -> profile "nndb")
    print(energy("GCACG", "CGUGC", "(((((&)))))", "nndb"), "expect -6.04")
    print(energy("AGCGCU", "AGCGCU", "((((((&))))))", "nndb", symmetry=True), "expect -7.94")
    print(energy("GGUCGUGU", "GCGUGGUC", "((((((((&))))))))", "nndb"), "expect -3.62")
