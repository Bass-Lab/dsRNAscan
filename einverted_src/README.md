# dsRNAscan einverted (self-contained fork)

`dsrnascan_einverted.c` + `emboss_minimal.h` are a **self-contained fork** of EMBOSS
6.6.0 `einverted` - not a patched copy. All EMBOSS library dependencies
(libajax/libnucleus/plplot) have been replaced by the `emboss_minimal.h` shim, so the
program builds with a plain C compiler and no EMBOSS install.

Because it is a fork (not stock-plus-a-patch), **`dsrnascan_einverted.c` is the single
source of truth.** The filenames were deliberately renamed away from upstream
(`einverted.c` -> `dsrnascan_einverted.c`) so it is clear at a glance that these are
modified files, not stock EMBOSS.

## What was changed from upstream

Grep the source for `MODIFIED` to see every deviation with a runtime effect:

1. **Standalone conversion** - EMBOSS calls satisfied by `emboss_minimal.h`. Pervasive
   and not individually tagged; it is why a clean patch against upstream is no longer
   possible.
2. **G-U wobble pairing** - tagged `MODIFIED (GU wobble)`.
3. **DP matrix reuse + traceback leak fix (v0.5.5)** - tagged `MODIFIED (reuse)`. The DP
   tables (`matrix`, `back`, `localMax`) are allocated **once** and reused across all
   input records; per record only the single wrap-around row (read stale at `i=0`) is
   cleared, instead of re-zeroing the whole O(maxrepeat^2) matrix per window. The
   alignment-traceback scratch is now also freed on the early-exit path (a leak that
   only mattered once dsRNAscan began batching many windows through one process). Output
   is byte-identical to the pre-reuse build.

`gu_wobble.patch` records **only** change (2) as a diff against the original upstream
`emboss/einverted.c`, kept for provenance. It is **not** applied by the build and does
**not** capture changes (1) or (3).

## Build

    ./build.sh                 # -> ./einverted for the current platform

Then copy the result to `../dsrnascan/tools/platform_binaries/einverted_<os>_<arch>`
(the shipped binary names are unchanged and are what setup.py / the runtime look for).
Build on each target platform, or cross-compile, to refresh all shipped binaries.

> **Release gate:** because batching runs many windows through one einverted process,
> the traceback leak fix in change (3) is mandatory on every platform - an unfixed
> binary will OOM on large scans. All platform binaries must be rebuilt from this
> source before a release, not just some.
