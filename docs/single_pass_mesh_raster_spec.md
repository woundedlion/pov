# Single-pass mesh rasterizer — design

Status: **REJECTED, measured 2026-07-18.** Do not build this. The go/no-go
probe (branch `probe/single-pass`, `cbeb26e4`) returned
`N_scheme/N_today = 0.83-0.96` against the 0.75 threshold this document set —
no solid clears it. The refutation is in "Why it fails" at the end; read that
before proposing any variant. The design and its two constraints are kept
because they are correct and they bound what any successor may attempt.

## The problem

`Scan::Mesh::draw` rasterizes each face independently through
`rasterize_face`. On IslamicStars' heaviest solids that spends **51 % of render
on SDF evaluations of pixels that are then rejected** (`d >= pixel_width`):
43.8 ms of an 85.5 ms frame, against 14.4 ms actually shading. Measured
41,820 probes/frame produce 18,604 shade events — a 2.25× probe-to-paint ratio
on top of a 1.79 shade-events-per-pixel overdraw.

## Two constraints, both learned the hard way

### C1 — adjacent faces disagree about their shared edge

`Face::distance` is a **gnomonic** distance in each face's own tangent plane.
For a point near the edge between faces A and B, `d_A` and `d_B` are computed
in *different* planes, so `d_A != -d_B`. On stretched faces the disagreement
exceeds the AA band.

This killed the per-pixel face-claim mask (`8cb9809b`, reverted `e7264a2c`):
letting face A claim a pixel it believes it fully covers dropped B's
anti-aliased contribution and painted the wrong colour, visibly roughening
wavefront edges.

**Therefore: no ownership. No claiming. No "each pixel drawn once."** A pixel
in a shared-edge band must be evaluated by *both* adjacent faces and blended,
exactly as today. The current look *is* that blend, disagreement included.

### C2 — per-face spans cannot beat a per-face bounding rectangle

The per-row span lever (`af49fdb6`, reverted `f7524d19`) computed each face's
true silhouette per row. It lost on device: 31.7 ms of span construction
against a 7.2 ms probe saving. The structural reason survives any tuning:
these faces are **13–20 columns wide**, a row's lit silhouette is **~5
columns**, and every emitted span must carry the AA fringe (1.25 col/side)
plus outward rounding — **~4.5 columns of overhead per span, often two spans
per row**. A per-face silhouette plus its mandatory fringe is simply not
smaller than the face's own bounding rectangle.

**Therefore: the fringe overhead must be amortized across the row, not paid
per face per row.** That is the one thing a whole-mesh pass changes.

## The design

Forward scanline over the whole mesh — *not* inverse per-pixel point location.

### Invariant that makes it safe

> The pass evaluates a superset of the (face, pixel) pairs that the current
> rasterizer evaluates **and does not reject**.

Skipped work is only ever an evaluation whose result would have been discarded
by `d >= pixel_width`. Therefore the framebuffer is **bit-identical** to
today's, and AA quality is not a judgement call — it is a diff. This is what
distinguishes the design from the claim mask, which changed output by
suppressing contributions.

### Stages

**S1 — static topology, once per shape spawn.** Mesh topology is fixed between
spawns; only the transform changes per frame. Build the unique-edge list with
each edge's adjacent face(s). `MeshState` carries flat `faces`/`face_offsets`/
`face_counts`, so this is a one-time O(I) pass, not per-frame work.

**S2 — per frame.** Vertices are already projected (`is_mesh_transform`). For
each visible edge compute its row range and one exact anchor crossing.

**S3 — per row, incremental.** Maintain an active-edge list. Each active
edge's azimuth crossing satisfies `cos(theta - psi) = k * cot(phi_y)`, with
`psi` and `k` row-independent. Rather than an `acos` per edge per row (the
~933 device cycles that sank the per-row lever), step it:

    du = k * (cot(phi_{y+1}) - cot(phi_y))      // cot from a row LUT
    dtheta ~= -du / sqrt(1 - u^2)

and **widen the span by the accumulated step-error bound**, re-anchoring with
an exact `acos` when the bound exceeds ~0.5 px. Conservative widening is legal
— the invariant needs a superset, not exactness. Target: 10–20 cycles per
edge-row against ~933 today.

**S4 — per row, assemble.** From the row's crossings, derive each active
face's span(s), padded by the AA fringe. Because edges are shared, the pad at
each crossing serves both adjacent faces — this is where C2's overhead
amortizes.

**S5 — per pixel.** For each face whose padded span covers the pixel, call
`Face::distance` and blend **exactly as `rasterize_face` does today**. One
evaluation in a face interior, two in a shared-edge band, more where three or
more faces meet at a vertex. Unchanged shader, unchanged alpha, unchanged
`filter_blend`.

### Expected yield

Per quadrant row (72 columns): ~131 face-pixel paints (the measured 1.8
shade-events/px, which is the AA itself and is irreducible) against ~290
probes today. The scheme's probes approach paints plus a per-crossing margin,
~140/row — **roughly −50 % probes, ~23 ms** on the worst shape. With the
landed clip cull that puts solid 11 near 61 ms against the 62.5 ms window and
the other five red solids comfortably green.

## Risks, in order

1. **Non-manifold and self-intersecting faces.** Registry solid 8 is 100 %
   self-intersecting (crossed hexagons), and hankin straps need not tile the
   sphere. The row is therefore **not** a clean partition. S4 must derive each
   face's spans from its own edges' crossings rather than assuming intervals
   alternate between exactly two faces; the sharing benefit degrades
   gracefully, correctness does not depend on it.
2. **ITCM.** The phantasm image has ~7 KB of RAM1 code headroom; the per-row
   span attempt ate 6.5 KB of it. This design is only viable if it **replaces**
   the per-face path for meshes rather than adding beside it, keeping code size
   roughly neutral. Verify early, not at the end.
3. **Construction cost.** The predecessor died here. S3's incremental stepping
   is the whole bet; if per-edge-row cost cannot be driven near ~15 cycles the
   design fails for the same reason.
4. **Vertex neighbourhoods.** Where 3+ faces meet, all incident faces must be
   evaluated in the fringe or AA degrades at every vertex — the most likely
   place for a subtle quality regression.

## Go/no-go before implementing

Do not build this before an offline probe answers, on the real solids:

- **Ceiling.** Simulate the scheme's evaluation set against the current one and
  report the ratio. If it is not near the predicted −50 %, stop.
- **Superset check.** Confirm the simulated set contains every (face, pixel)
  pair the current renderer evaluates and does not reject. This is the
  bit-identical claim, verified before any rasterizer is written.
- **Construction budget.** Price S3 per edge-row and multiply out. If it is not
  comfortably under the saving, stop — that is exactly the measurement the
  per-row lever skipped.

Predecessors and their causes of death: `e7264a2c` (claim mask — AA, C1),
`39db0225` (cap-span clip — no-op, and its landing figure came from comparing
`perf_bench` across sessions), `f7524d19` (per-row spans — construction cost,
C2).


## CORRECTION (audited 2026-07-18)

The verdict below stands — `N_scheme/N_today = 0.83-0.96` reproduces exactly —
but **both of its supporting claims are wrong.** Audit branch `hs-wt-audit`
off `cbeb26e4`.

**The "91-98 % overlap, irreducible under C1" claim is inflated 3-4x.** True
C1-protected fringe-vs-fringe overlap is **22-31 % of `N_today`**; the rest is
*bounding-rectangle* overlap on rows that declined to the extent (64-86 % of
the scheme's evaluations sit on such rows, overlapping at 83-96 %). Rectangle
overlap is not protected by C1 — it is waste. Cross-check without any scheme:
1.48-1.73 contributing faces per painted pixel, 54-70 % of contributing pairs
on a multi-contributor pixel.

**Rows do not decline because of ~4.5 columns of overhead, nor because of
multi-span star rows.** Decline rate conditioned on span count is flat (solid
11: 79.7 % at one span, 76.2 % at two), and mean spans/row is 1.83 on emitted
vs 1.84 on declined rows — no signal. The cause is that the AA fringe converts
to azimuth as

    reach_theta = FRINGE_PX * (2*PI/W) / sin(phi)

and **the rendered band is a polar cap** (phi ~ 0-45 deg). Measured pad per
side: **18.98 columns at row 3**, 4.80 at row 12, 2.48 at row 24, 1.76 at row
36. The 1.25 col/side figure used throughout this document holds only at the
equator. Decline rate tracks it monotonically (90.6 % at row 3 to 66.5 % at
row 36); rows y<=2 cannot build spans at all. Declined-row spans measure
75-85 % fringe stadia, not interior. A secondary effect compounds it: declines
concentrate on faces whose extent is already tight (13.74 declined vs 19.00
emitted on solid 11). Not a quadrant artifact — full-frame shows the same
mechanism, milder.

**Consequence — a latent correctness bug in the CURRENT renderer.** The
whole-face extent path pads by a constant `1.25 * (2*PI/W)` regardless of
latitude (`Face::get_horizontal_intervals`), when covering a one-pixel arc
fringe requires `1.25/sin(phi)` columns. Master therefore **clips AA fringe
near the poles**. This is almost certainly the unexplained observation that
per-row spans painted 1-19 more pixels/frame than the extent path. Independent
of any performance work.

**What the audit did NOT overturn:** the floor `N_ideal/N_today = 0.41-0.46`
is real, and 27-43 % of probed face-rows contribute nothing (26-35 % of
`N_today`). The best legal variant found — intersecting row spans with the
whole-face extent, which can only shrink a row and never declines — reaches
**0.78-0.92** with zero misses, still short of the 0.75 threshold and without
pricing construction. Rows that cannot produce spans carry 23-29 % of
`N_today` on every solid. One corner remains unmeasured: the stadium's
arc-clip extension uses a worst-case orientation bound
`pad = REACH * max(1/sin phi_row, 1/sin phi_arc)`; tightening it is estimated
(not measured) to move solid 11's declined-row span from ~29.4 to ~22 columns,
which would still not flip most declines.

## Why it fails — as originally argued (superseded above)

Measured over solids 11/6/8/19/9/7, 64 poses, quadrant clip, with the probe's
replica of today's run construction reproducing `pixels_tested` exactly.

**The geometry does admit the predicted win, and the fringe eats it.** With a
*zero* fringe pad the minimal per-(face,row) arc containing that face's
contributing columns sits at `N_ideal/N_today = 0.43-0.46` — essentially on
the contributing floor (`N_contrib/N_today = 0.41-0.46`). Adding the mandatory
AA pad and outward rounding takes it to 0.83-0.96. **The fringe consumes about
three quarters of the available win.** That is C2 quantified, and it is not a
tuning failure: 47-80 % of rows have their padded spans decline to the
whole-face extent because the padded span is not narrower than today's
rectangle.

**The whole-mesh pass cannot fix that, and C1 is the reason.** 90.9-98.2 % of
the scheme's evaluations land on pixels covered by two or more faces' padded
spans — the fringe overlap is nearly total. But C1 forbids collapsing it: each
adjacent face must still evaluate its own padded span, because the AA blend
*is* their disagreement. So sharing an edge between two faces amortizes
**construction** (a real 1.58-1.60x reduction in edge-rows) and **never the
pair count** — and the pair count is the only thing the probe saving depends
on. `N_scheme` is therefore identical in kind to the per-row lever's
evaluation set, which `f7524d19` already measured and rejected.

**The central premise of this document — that a whole-mesh pass changes the
probe count — is false by C1's own logic.**

Independently, S3's incremental stepper does not work either. It measures
2.37-2.83 ns/edge-row against 2.14-2.44 ns for an exact `fast_acos`: a 1.11-
1.24x **regression**, because the error-bound arithmetic costs more than the
transcendental it avoids, and it still re-anchors 24-29 % of the time at a
correctly-calibrated 0.5 px bound. It also attacks the wrong 3 %: `acos` is
2.3 ns of the 79-89 ns/edge-row that `emit_row_spans` costs. The predecessor's
~933 cycles is band/clip arithmetic, crossing pairing, insertion sort and the
depth walk — the ~15 cycle target was never reachable.

## What the probe vindicated

- **The superset invariant HOLDS.** Zero contributing pairs missed, pair-by-
  pair, every solid x 64 poses. The bit-identical claim was sound; there was
  simply nothing worth buying with it.
- **Vertex neighbourhoods are clean.** 68k-227k contributing pairs lie within
  1.5 px of a projected vertex (~40 % of all contributing pairs on solid 11 —
  these meshes are vertex-dense) and **zero** are missed. The predicted AA-loss
  risk at 3+ face junctions does not materialize.
- **Self-intersecting faces are clean.** Solid 8 (100 % crossed hexagons) holds
  the superset exactly and sits mid-pack at 0.861, so non-manifold geometry was
  never the obstacle.

## Where to look instead

Not span tightening. Four attempts have now died on the same structural fact,
and C1+C2 close the family: no per-face-span scheme can reduce the pair count
below the fringe, and the fringe is the AA. What remains is **coverage
reduction** (fewer or larger faces — an art decision) or **per-probe cost**.

On per-probe cost, one figure deserves a look before it is assumed harvested:
the device spends ~1,878 scan cycles per shaded pixel at ~2.25 probes per
shaded pixel, implying **~830 cycles per probe**. That is implausibly high for
a sector walk (binary search plus three edge distances) and suggests the
per-pixel loop carries substantial cost outside `Face::distance` itself. Note
also that no probe count has ever been measured on device — `HS_SCAN_METRICS`
is not compiled into the profile image — so every probe figure in this
document's lineage is host-simulated and the absolute baseline (41,820
probes/frame on solid 11) does not reproduce at quadrant clip (the probe
measures 19,330). Reconcile that before reusing it.
