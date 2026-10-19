# Turner 2004 duplex energy evaluation (reference for the Rust rescore)

`turner2004.json` (NNDB tables, kcal/mol; layout in `meta.layout`), `duplex_energy.py` (`energy(s1, s2, structure)`),
`validate.py` (writes `validation_report.txt`, `validation_results.tsv`), `build_turner2004.py`, `nndb_raw/` + `SOURCES.md`.
Evaluation only, no folding DP. ViennaRNA source was read for conventions; no code or parameter file was copied.

## What `RNA.duplexfold`'s energy contains  (V = tested, S = read in duplex.c, consistent with data)
- E = init + two helix-end terms + stack/bulge/interior loops between consecutive pairs (V).
- End term = terminal mismatch if both flanking nucleotides exist, else the 3' or 5' dangle that exists, + AU/GU end penalty (V).
  Flanks are neighbours in the full sequences, so the returned structure keeps one `.` per side when available. Evaluating
  without them changes 80% of energies (mean +0.97, max 3.2 kcal/mol).
- `fold_compound(t1&t2).eval_structure(structure without '&')` on the trimmed strands *with flanks* equals duplexfold's energy
  for 2000/2000 pairs. `slice_duplex()` makes that trim from `result.i/.j`.
- Parameters are the 0.1-rounded NNDB values: stacks, init 4.1, AU/GU end 0.5. 1-nt bulge = 3.8 + stack (no end penalty);
  larger bulge = init + 0.5 per AU/GU closing pair; interior loop: 1x1/1x2/2x2 tables, else init + min(3.0, 0.6|n1-n2|) +
  first mismatches (none for 1xn) + 0.7 per AU/GU closing pair; <=30 unpaired per loop (S); lonely pairs allowed.
- N is scored as the worst of A/C/G/U for its term. 37 C only. `dangles` 1/2/3 give identical energies; 0 differs (not modelled).

## Profiles
`duplexfold` (default) reproduces the above. `nndb` follows the NNDB text: 2-decimal stacks, 4.09, 0.45, special C bulge,
bulge state count, GGUC/CUGG motif, uncapped asymmetry; the NNDB worked examples reproduce (-6.04, -7.94, -3.62).
Values in `duplexfold` that are not in NNDB were measured against `eval_structure`, not copied from a file: asymmetry cap 3.0;
generic first-mismatch GG -1.0, UU -0.6 (NNDB -1.2, -0.7); 2x3 GG -0.7, UU -0.3 (NNDB -0.8, -0.4); 2x3 GA with AU/GU closure -1.1
(NNDB -1.2). Whether measured numbers count as derived from ViennaRNA's parameters is a licensing call; `nndb` avoids it but does not match.

## Agreement (2000 pairs of 15,366 chr19 candidates: 20-60 nt n=748, 61-499 n=956, 500+ n=296; duplexes to 1386 bp)
**Table 1. Agreement of the evaluator with `RNA.duplexfold` energies (2000 chr19 candidate duplexes).**

| profile vs `RNA.duplexfold` | <=0.01 | <=0.1 | <=0.5 | mean abs | max abs |
|---|---|---|---|---|---|
| `duplexfold` | 100% | 100% | 100% | 0.000 | 0.00 |
| `nndb` | 1.5% | 8.5% | 47.9% | 2.05 (500+ nt: 7.87) | 69.0 |

*Methods.* 15,366 arm pairs from `dsrnascan-core` on `validation/out/chr19_300kb.fa` (`--no-select --strand both --min-match-perc 70`); 2000 sampled
(seed 1) by `len(i_seq)`: 20-60, 61-150, 151-499 and 500+ nt strata, the 500+ pool (296) taken whole and its shortfall added to 20-60. Per pair,
`RNA.duplexfold(i_seq, j_seq)` (ViennaRNA 2.7.0, 37 C, dangles 2), strands cut to the structure including flanks, `energy()` under each profile;
diff = mine - duplexfold in kcal/mol, rounded to 0.01. Columns: fraction of pairs with |diff| <= tolerance, mean |diff|, max |diff|. Per-stratum rows
(all 100% for `duplexfold`) are in `validation_report.txt`. Script `validate.py`; per-pair output `validation_results.tsv`.

`duplexfold` is also exact on ~84,000 constructed duplexes (every 1x1/1x2/2x1/2x2 entry, every end pair x flank, generic/2x3/1xn mismatch
rules, bulges 1-30, large loops, N). Its 10 worst pairs have diff 0.00 (no responsible term). The 10 worst for `nndb` are 500+ nt arms with a
20-29 nt asymmetric loop (e.g. 3x26: 16.7 vs 5.9); term = uncapped asymmetry penalty.

## NNDB text -> duplexfold, mean effect over all pairs / over affected pairs (kcal/mol)
stacks 2-dec -> 0.1: -0.96 / -0.96; asymmetry cap: -1.47 / -6.9 (21% of pairs); interior mismatch values: +0.33 / +0.64; no bulge state
count: +0.41 / +1.0; end penalty 0.45 -> 0.5: +0.21 / +0.23; no special C bulge: +0.12 / +1.0; no GGUC motif: +0.09 / +4.9; init: +0.01.
The self-complementary symmetry term (+0.43) is available as `symmetry=True`; duplexfold never applies it.

## Parity caveats
Structures are not reproduced: duplexfold picks the end pair by strict `<` (i ascending, j descending) and backtracks to the first matching
predecessor (k descending, l ascending), so ties extend the helix (S, untested). `percent_paired` and helix length depend on this.
Use integer centi-kcal in Rust. Loops >30 nt never occur in duplexfold output (log extrapolation checked to 0.01 only).
