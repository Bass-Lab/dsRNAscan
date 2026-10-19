# dsrnascan-core design

One job: find intramolecular dsRNA candidates (inverted repeats within `max_span`) in genomic
sequence, score them with nearest-neighbor energies, and hand them to Python for ML scoring and
output. One sweep of the sequence; everything else reads what the sweep recorded.

## 0. Default: dsRNAscan 0.5.x's output, from one sweep

The default reproduces 0.5.x's `-w 10000 -s 500` output **exactly**. On chr19:1,000,001-1,300,000,
both strands, the result is identical to 0.5.x in every column except the corrected
`longest_helix`:
- 725/725 dsRNAs at `-s 500`, and 840/840 at `-s 250`.
- Zero differing cells (coordinates, score, match counts, gaps, dG, base pairs, percent_paired,
  trimmed coordinates, sequences, structure).

How:
- **Windows do not change DP values.** A cell depends only on cells with a larger 5' position
  and a smaller 3' position, which lie in the same window. So every cell einverted computes in a
  window has the value the single sweep computes.
- **Windows change only einverted's reporting.** It keeps one best cell per row among partners
  inside the window, then runs its slot bookkeeping (`back`, `localMax`, `lastReported`).
- **The sweep records each row's running maximum** (`scan::RowSteps`), so the row's best within
  any window is a lookup.
- **einverted's report loop is replayed per window on those rows** (`einverted.rs`, ported from
  EMBOSS einverted with dsRNAscan's G-U patch). One edge case is matched: einverted never reports the
  adjacent-base cell (column 0). It reports the next cell and keeps that match in the score
  (`COL0_META`).
- **The 0.5.x pipeline then runs on the union:** match pre-filter, fold (`duplex.rs`),
  percent_paired / min_bp filters, then 0.5.x's nested elimination (`select::nested_05x`).
- **Lock-in test** (`validation/single_window.py`): on sequences exactly one window long the replay
  reports exactly einverted's hits. 80/80 windows, 225/225 hits; also 10/10 at `--gap 4`.

Speed on that region (10 threads; 0.5.x `-c 8`: 39 s wall, ~252 CPU-s):
- exact fold: 9.8 s wall, 66 CPU-s.
- `--banded`: 3.9 s wall, 26 CPU-s. Same 725 dsRNAs; fold values differ in 2.

**Parallelism:** one rayon pool is shared by sequences x strands, sweep blocks, windows and folds.
Output is byte-identical at any thread count.

Sections 3 and 5 describe `--windowless`: selection straight from the sweep, with no windows. It
finds arm-sharing partners that 0.5.x cannot. Its volume is controlled by greed / rescue /
`--compact`.

## 1. Search space

Every candidate base pair is a cell in pair space: row `q` = position on the 3' arm, column `k`
pairs it with 5'-arm position `p = q - 1 - k`. A duplex is a path that grows outward from its
innermost pair (the origin) to its outermost pair, so it covers a box in the dot plot:
5' arm `[p_out, p_in]` x 3' arm `[q_in, q_out]`.

Only `k < max_span - 2` is computed (einverted's `-maxrepeat` semantics), so the work is
`sequence length x max_span` cells, each computed exactly once.

## 2. Scoring (v0: einverted-compatible)

Smith-Waterman with a linear gap, the recurrence and tie rules of einverted (EMBOSS 6.6.0, G-U
counted as a match): match +3, mismatch -4, gap 12, threshold 75. Ties prefer a 5' bulge, then
a 3' bulge, then extending the helix. On the same sequence the DP values are identical to
einverted's.

Planned: replace the match/mismatch table with a nearest-neighbor stacking table, so GU-rich or
AU-rich repeats are not scored as if they were GC helices. This is a table swap, not a redesign.

## 3. What counts as a candidate (`--windowless`)

Every cell carries its path's origin and gap count (packed into one u64). For each origin we keep
its best cell: highest score, then the later row, then the first column, which is einverted's
preference. The match percentage is a constraint on that choice, not a filter applied afterwards:
the best cell is the best one whose path meets `min_match_perc`. Otherwise a path that extends
into low-identity flanks would fail the cutoff and take its high-identity core down with it.

Origin, gap count, score and box together give matches, mismatches and gaps exactly. The forward
tie rules equal einverted's traceback rules, so no traceback pass is needed.

einverted keeps one partner per 3' position and one 3' position per 5' start. That is why it drops
a weaker hit that shares an arm with a stronger one, and why dsRNAscan <= 0.5 needed overlapping
windows. Per-origin bookkeeping removes the cause.

## 4. Parallelism without duplicates

Each sequence is cut into blocks (automatic: about two per thread, 4 halos to 4 Mb long). A block
sweeps `max_span + 2` rows before (and, for `--windowless`, after) the rows it owns, and reports
only what belongs to rows it owns.
- **Exact:** a path's rows all lie within `max_span` of its origin, so every owned path is
  computed exactly as in a single sweep.
- **Each candidate reported once:** each origin is owned by one block.
- **Verified:** block-parallel output is byte-identical to a single sweep.

## 5. Selection (`--windowless`)

Hits are kept strongest-first. A hit is dropped when at least `overlap_frac` of BOTH of its arms
lies inside the corresponding arms of a hit already kept.
- Default `overlap_frac = 1.0` drops only fully nested hits, matching 0.5.x's nested elimination.
- `--compact` (`0.5`) also drops hits that substantially overlap a stronger one on both arms.
  Smaller output, and fewer folds.

**Greed (default 1.0).** A hit is also dropped when one of its arms overlaps the SAME arm (5'/5' or
3'/3') of a stronger kept hit by at least half its length. Under `--greed G` it is dropped only if
it scores below G x the stronger hit's score.
- Cross-role sharing is never suppressed, as in einverted. That is an arm used as one hit's 5' arm
  and another's 3' arm, e.g. alternating Alus.
- `--greed 0` keeps every arm-sharing hit. On chr19 that gave 13x 0.5.x's output, mostly
  near-duplicates.

**Rescue (default on; `--no-rescue` turns it off).** This reproduces what 0.5.x's overlapping
windows did. A weaker partner of a shared arm survives when the stronger hit's other arm lies
outside the weaker hit's span, because some window held the weaker pair without the stronger one.

chr19 300 kb, against 0.5.x (725 dsRNAs, 353 scored likely edited by the stability model):

| mode | dsRNAs | 0.5.x recalled | 0.5.x likely-edited kept | new likely-edited | window* |
|---|--:|--:|--:|--:|--:|
| greed 1 + rescue (default) | 2,207 | 612 | 285 | 328 | 53 |
| greed 1 + rescue + compact | 2,027 | 605 | 283 | 311 | 41 |
| greed 1, no rescue | 1,141 | 353 | 150 | 120 | 27 |
| greed 0 | 9,691 | 700 | 337 | 2,061 | 364 |

\* dsRNAs touching chr19:1,028,736-1,030,616 (0.5.x: 12). None of these rule sets reproduces 0.5.x
exactly: 0.5.x's set depends on which partner einverted ranked best inside each window.

Ordering: the quality filters (match %, then the energy-derived ones) run BEFORE selection, so a
hit that fails a filter can never suppress one that passes.

## 6. Strands and orientation

- **Scanning.** The minus strand is scanned as the complement of the forward sequence, not
  reversed, so coordinates stay on the forward axis. Pairing on the complement is exactly
  minus-strand RNA pairing, including wobble: a minus-strand G-U sits on genomic C/A.
- **Coordinates.** 1-based inclusive on the forward axis. `i` = lower-coordinate arm,
  `j` = upper (BED-friendly, unchanged from 0.5.x).
- **Sequences and structure.** Always 5'->3' along the transcribed strand. On `-` the arms are
  reverse complements, and the transcript's 5' arm is `j`, not `i`.
  - 0.5.x folds and reports `i&j` on both strands, which puts the 3' arm first on `-`.
    Energetically equivalent, but backwards for visualisation.
  - Proposed: report `arm5`/`arm3` sequences and a structure in transcript order on both
    strands. **Decision pending.**
- **Duplex-trimmed coordinates.** On `-`, `eff_*` coordinates map trim positions from the reversed
  arms back to the forward axis (as 0.5.x does: `eff_i_start = i_end - i_trim_end + 1`, etc.).

## 7. Output contract (draft)

Required, per the dsRNAscan manuscript and ML models:
- `chrom, strand, i_start, i_end, j_start, j_end`: search extent
- `eff_i_start, eff_i_end, eff_j_start, eff_j_end`: duplex-trimmed extent
- `dG`: duplex free energy, kcal/mol
- `percent_paired`: 2 x pairs / trimmed duplex length x 100
- `longest_helix`: longest run of directly stacked pairs. This is what the training data's
  `bp_count` is: 99,999/99,999 joined chr1 dsRNAs match. The 0.5.x `find_longest_helix` does
  not compute this; see the validation notes.
- `length`: shorter arm, from search coordinates. `loop_length`: `j_start - i_end`.
- Arm sequences and dot-bracket structure, transcript order (see 6).
- ML scores.

Optional (debug): `score`, `match_perc`, `gaps`.

## 8. Energy rescore

`src/duplex.rs` is a port of ViennaRNA's `duplexfold` (duplex.c, 2.7.0; credit Lorenz et al. 2011).
- **Same as ViennaRNA:** recursion, integer arithmetic, tie rules and flank handling. Energy,
  structure and `.i`/`.j` are identical to `RNA.duplexfold` on 3500/3500 test pairs.
- **Helix and pair counts:** `longest_helix` and the pair count come from the traceback (a stacked
  step extends the helix, a loop resets it), so no dot-bracket parsing.
- **Band option:** `duplexfold_band` folds only an anti-diagonal band around the search path.
  - The path itself lies within `gaps` of its end diagonals, so `path_band` adds `margin` to that.
  - margin 30: 99.25% identical, 4x faster. margin 60: 99.75% identical, 2.7x faster.
  - The rest are hits where duplexfold's best duplex is in another register: tandem repeats, and
    Alu pairs shifted by about one monomer.
  - **Default is exact** (matches published output). `--banded` switches the band on.
- **Exact pruning:** the full fold skips any loop whose lowest possible energy cannot beat the
  current best. Results are unchanged (3500/3500 identical) and it runs 2x faster than ViennaRNA.

Cutting rescore volume without changing results: rescore lazily inside selection. Walk candidates
strongest-first. Skip any already covered by a kept hit (no rescore needed). Rescore the rest, and
keep a hit only if it passes the energy filters. This is equivalent to filter-then-select, but
folds far fewer candidates.

Parameters and a reference evaluator are in `energy/` (see its README).
- **Source:** Turner 2004 values taken from NNDB, not from ViennaRNA's files, because its no-fee
  clause is likely GPL-incompatible.
- **Profile `duplexfold`:** NNDB values plus six conventions measured against ViennaRNA's output
  (asymmetry cap 3.0; four first-mismatch values; 0.1 rounding). It reproduces `RNA.duplexfold`
  energies exactly: 2000/2000 chr19 pairs, 535/535 on an independent chr1 region.
- **Profile `nndb`:** NNDB taken literally. Mean difference 2.05 kcal/mol, up to 69 for long arms
  with asymmetric loops.
- **Flanks:** the evaluator must be given the flanking nucleotide on each side.
- **Licensing decision pending:** whether the measured conventions are acceptable.

## 9. Validation (validation/, run from that directory)

| script | checks |
|---|---|
| `synthetic.py` | one arm with two partners, every layout; both pairs must be found |
| `parity.py <fa>` | every einverted single-pass hit appears among core candidates with identical coords, score and counts; block-parallel output is identical |
| `vs_pipeline.py <fa> <old merged_results>` | end-to-end against dsRNAscan 0.5.x (same RNAduplex filters, then the box rule) |

## 10. Later

PyO3 + maturin wheels (Python keeps the CLI, output formats and ML); fusing both strands into one
pass (they differ only at wobble pairs); SIMD kernel; streaming selection so candidates never
accumulate per chromosome.
