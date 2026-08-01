# Design note: the SDF shape interval protocol has no per-scan setup hook

Status: **recorded, not built.** A design change to `core/render/sdf.h`, sized
and scoped here so a future session can decide against real numbers instead of
re-deriving the shape of the problem.

## The protocol today

`scan_region<W, H>` (`core/render/scan.h`) drives a shape through exactly one
entry point per row:

```cpp
bool get_horizontal_intervals<W, H>(int y, OutputIt out) const
```

There is no call the shape sees *once per scan*. Every per-scan quantity must
therefore be recomputed per row, hoisted into the constructor (which does not
know `W` or `H`), or memoized behind mutable state.

## What that costs per row

Two classes of work repeat on every row of every scan:

- **Scratch setup in the CSG composites.** `Union`, `SmoothUnion`, `Subtract`
  and `Intersection` each open a `ScratchScope` on `scratch_arena_b` and
  placement-new their span buffers on entry, then tear them down on return.
  Composites nest, so a two-level tree pays this at each level, per row.
- **LUT initialization re-tests.** Six sites re-test
  `TrigLUT<W, H>::initialized` per row, although `init_geometry_luts<W, H>()`
  runs from engine setup before the first frame and all three drivers go through
  it. The guards remain only as a lazy fallback for unit tests and offline
  tools.

Neither is expensive on its own; both are pure per-row overhead that a per-scan
hook would pay once.

## What a hook would look like

Add an optional `begin_scan<W, H>()` that `scan_region` calls before the row
loop and that composites forward to their children, with a `requires` probe so
shapes that need nothing stay unchanged. The scan-lifetime scratch scope then
lives in `scan_region` alongside its own `intervals`/`norm` buffers, and the
`TrigLUT` guard collapses to the one call already made there.

## Why it is not built here

The arena is a LIFO bump allocator shared with `Pixel::Feedback::flush`, and the
existing per-row scopes are what keep a composite's span buffers from outliving
the row that produced them (`core/render/sdf.h`, the scratch-arena contract
comment). Lifting them to scan lifetime changes when every CSG node's scratch is
reclaimed, across a nested tree — a correctness-sensitive refactor of the
deepest render chain, not a patch. It wants its own measurement pass on device
(`docs/profiles`) to establish that the per-row setup is worth the change at
all: `scan_face_setup` measurements in `docs/probe_path_open_items.md` are a
reminder that static instruction counts predict this path poorly.
