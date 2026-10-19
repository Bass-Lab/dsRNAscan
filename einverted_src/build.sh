#!/usr/bin/env bash
# Build the self-contained einverted used by dsRNAscan.
#
# This is a de-EMBOSS'd standalone (no libajax/libnucleus/plplot needed) with the
# G-U wobble patch already applied, plus the v0.5.5 "matrix-reuse" change: the DP
# tables are allocated once and reused across input records (only the single stale
# row is cleared per record) so that scanning many windows in ONE process no longer
# re-zeros the O(maxrepeat^2) matrix per window. Output is byte-identical to the
# previous einverted build.
#
# Produces one binary per platform; run on each target (or via CI) and copy into
# dsrnascan/tools/platform_binaries/.
set -e
cd "$(dirname "$0")"
OUT="${1:-einverted}"
gcc -O2 -o "$OUT" dsrnascan_einverted.c -lm
echo "built $OUT"
