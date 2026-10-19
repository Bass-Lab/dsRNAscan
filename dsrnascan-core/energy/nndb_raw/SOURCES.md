# NNDB sources (RNA, Turner 2004 set)

Downloaded 2026-10-17 (UTC) with curl from the Nearest Neighbor Database, https://rna.urmc.rochester.edu/NNDB/
(Turner & Mathews 2010 NAR 38:D280; parameters: Mathews et al. 2004 PNAS 101:7287, Xia et al. 1998
Biochemistry 37:14719 for Watson-Crick stacks). Files are stored unmodified; `../build_turner2004.py` parses them
into `../turner2004.json`. Base URL below: `B = https://rna.urmc.rochester.edu/NNDB`.

## Parameter tables (`params/`, used by the parser)

| File | URL | Table | JSON key |
|---|---|---|---|
| `turner_2004_watson_crick_stack.htm` | B/parameter_tables/turner_2004_rnastructure/turner_2004_watson_crick_stack.htm | WC nearest-neighbour stacks, dG37 (2 decimals), intermolecular initiation 4.09, per-AU-end 0.45, symmetry 0.43 | `stack`, `init` |
| `turner_2004_gu_pair.htm` | B/parameter_tables/turner_2004_rnastructure/turner_2004_gu_pair.htm | stacks containing GU, GGUC/CUGG special case, per-GU-end 0.45 | `stack`, `stack_special` |
| `turner_2004_watson_crick_stack_dg.txt` | B/parameter_tables/turner_2004_rnastructure/turner_2004_watson_crick_stack_dg.txt | all 36 WC+GU stacks, plain text, 0.1 kcal/mol (cross-check of the html values; used as is for `stack_1dp`) | `stack_1dp` |
| `turner_2004_terminal_mismatch_dg.txt` | B/parameter_tables/turner_2004_rnastructure/turner_2004_terminal_mismatch_dg.txt | terminal mismatch (exterior helix end with both flanks) | `terminal_mismatch` |
| `turner_2004_dangle_dg.txt` | B/parameter_tables/turner_2004_rnastructure/turner_2004_dangle_dg.txt | 3' and 5' dangling ends | `dangle3`, `dangle5` |
| `turner_2004_loop_dg.txt` | B/parameter_tables/turner_2004_rnastructure/turner_2004_loop_dg.txt | bulge and interior loop initiation by size 1-30 (hairpin column ignored) | `bulge.init`, `interior.init` |
| `turner_2004_bulge_loop.html` | B/parameter_tables/turner_2004_rnastructure/turner_2004_bulge_loop.html | bulge initiation (n=1 listed as 3.81), special C bulge -0.9 | `bulge` |
| `turner_2004_internal_loop_dg.txt`, `turner_2004_internal_loop.htm` | B/parameter_tables/turner_2004_rnastructure/turner_2004_internal_loop_dg.txt (and `.htm`) | interior loop asymmetry 0.6, AU/GU closure 0.7, 2x3 and generic first-mismatch bonuses | `interior` |
| `turner_2004_int11_dg.txt` | B/parameter_tables/turner_2004_rnastructure/turner_2004_int11_dg.txt | 1x1 interior loops | `int11` |
| `turner_2004_int21_dg.txt` | B/parameter_tables/turner_2004_rnastructure/turner_2004_int21_dg.txt | 1x2 interior loops | `int21` |
| `turner_2004_int22_dg.txt` | B/parameter_tables/turner_2004_rnastructure/turner_2004_int22_dg.txt | 2x2 interior loops | `int22` |

## Documentation pages (formulas and conventions; not parsed)

B/rna_2004/rna_2004_watson_crick_helices.html, rna_2004_gu_pairs.html, rna_2004_dangling_ends.html,
rna_2004_terminal_mismatches.html, rna_2004_bulge_loops.html, rna_2004_internal_loops.html,
rna_2004_exterior_loops.html.

Not needed for a two-strand duplex and not downloaded/kept: hairpin loops, multibranch loops, coaxial stacking.

## SHA-256

```
fde6b078657e7511cbe3127426323665aa75da034a13c9670930b719834475c9  rna_2004_bulge_loops.html
5a57d2003e07de1cb8bc4b73d67401be0c86a36277df74cb8649366b330566f6  rna_2004_dangling_ends.html
ae198acadf54c87531cf156393d66dd207fcb2ee88e44a98fa74972957cd3339  rna_2004_exterior_loops.html
220bc452222015d7feb127f2494aedb962e896a3be0634209f1cc117f88634d4  rna_2004_gu_pairs.html
217978c4f4f0a0a10fe950e3fcf701f7bb146a620b90edaccbe9b175750bbebc  rna_2004_internal_loops.html
093959ef9a1c91c5eb8582dd502619b0a6cb0384a1c8239350288f272e328b26  rna_2004_terminal_mismatches.html
ac8658a8d456af7e9c3049afb2238225ee7b6e31fda486a2ced95d289d7749b8  rna_2004_watson_crick_helices.html
10981c3d24ce0717afb3a95770c9a9613cc608dddfd52ef3f223153b96c90de2  params/turner_2004_bulge_loop.html
a91e533c756525230d14696fe26c0dd9c3a3099bdb91cb3ada22af6f4a43ef88  params/turner_2004_dangle_dg.txt
ed426af901d921dd0069e4face2e15aa42580aba7d2966c7260427e73d7b35a6  params/turner_2004_gu_pair.htm
ccbb8992fc83b7476b7ac851d3a5db5c01a0de2fea7eedd72dd2e71f95cf1c66  params/turner_2004_int11_dg.txt
39f89b322360be5e76ddc16e2c02d1fde8a2004e1d42110f4632318b4e47e4db  params/turner_2004_int21_dg.txt
56ac7b6eae811db7e70a00443a264edaedfd9d914c110bb803b4e6714e52de9d  params/turner_2004_int22_dg.txt
faa2baa056fa58cd413dd77ee69015a06b9b0f31142ac6d8339701619a588af7  params/turner_2004_internal_loop.htm
e098e3b37c6980fb239b6bf0450de067d315a155051f54a252ea3ac26e7d3674  params/turner_2004_internal_loop_dg.txt
c1a330a3382e2f2d3dc26ba067b110d5c6c2cc3aa2986d174bcabe602c2ce07e  params/turner_2004_loop_dg.txt
0acbfbb32890b4ad55ca4193e755dceb6083bdbae6c9bc67a5299bb5c069f5b9  params/turner_2004_terminal_mismatch_dg.txt
2b1a9dccf1b2c80b9d3fb0b11251468082337dde261913e087f87de6a59aa263  params/turner_2004_watson_crick_stack.htm
5ffe3e1b9bc8fbcddfdfbcb67df0738e7a935cd85d889a18b3de62337915231e  params/turner_2004_watson_crick_stack_dg.txt
```

## Notes

* The NNDB site was redesigned; the pages are now under `/NNDB/rna_2004/`, the old `/NNDB/turner04/` paths return 404, and
  the `turner_2004.zip` link on the download page returned 404 on 2026-10-17, so the individual files above were used.
* The 0.1 kcal/mol `*_dg.txt` tables are the 2-decimal published values rounded to one decimal (checked for all 36 stacks in
  the parser). Intermolecular initiation (4.09) and AU/GU end penalty (0.45) have no plain-text file in this set; the
  `*_1dp` scalars (4.1, 0.5) in `turner2004.json` are those numbers rounded half away from zero (a derivation, not an NNDB listing).
