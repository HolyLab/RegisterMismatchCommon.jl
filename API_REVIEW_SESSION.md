# Session Handoff — 2026-05-18

## Plan
API_REVIEW_PLAN.md — RegisterMismatchCommon, v1.0.1

## What was just completed
CHUNK-002: export-mismatch-apertures
Added `mismatch_apertures` to the `export` list at line 25 of `src/RegisterMismatchCommon.jl`. Verified `:mismatch_apertures in names(RegisterMismatchCommon)` returns `true` after restart; 0 new ambiguities.

## Key decisions / shim choices
- Appended to the third export line (line 25) alongside `register_translate`.
- No tests needed: export visibility is structural, not behavioral.

## State of the codebase
- Files modified: `src/RegisterMismatchCommon.jl` (line 25, added `mismatch_apertures` to exports)
- Test suite: not re-run this chunk (export-only change; baseline was 140/140)
- Ambiguity count: 0 (unchanged)
- Staged but uncommitted: yes (`src/RegisterMismatchCommon.jl`)

## Cluster status
- mutating-pairs: 0 of 2 complete
- thresh-migration: 0 of 2 complete

## Next chunk
CHUNK-003: mismatch-kwargs-forwarding — replace explicit `normalization=:intensity` keyword in both `mismatch` methods with `kwargs...` forwarding to the downstream protocol.

## Watch out for
- The `normalization` keyword moves from this package to the downstream. Callers who mis-spell `normalization` will now get an error from downstream rather than from this package. Acceptable per the stated values.
- Both `mismatch` method signatures (the `{T<:AbstractFloat}` one and the promoting one) need updating.
