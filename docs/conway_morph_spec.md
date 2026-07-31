# ConwayMorph — Animated Conway-Operator Transitions for HankinSolids — Design Spec

**SUPERSEDED.** This document records the design as it landed and no longer
tracks the tree: the symbol names, module layout and tuning constants below have
since moved on. The current design of record is `docs/opchain_morph_spec.md`;
the source of truth for behaviour is `core/animation/opleg.h` and
`core/mesh/hankin.h`.

Status: LANDED 2026-07-16 (`06c95e93..190a91d8`). Replaces the
nearest-vertex-slerp `Animation::MeshMorph` + stochastic dissolve between
HankinSolids solids with true geometric morphs: each transition sweeps the
continuous parameter of the destination solid's own Conway chain, so faces
visibly truncate, expand, and twist from one polyhedron into the next. Exactly
one mesh is on screen at all times — there is no crossfade or dissolve
anywhere in the design (§2.5). `Animation::MeshMorph` itself remains for other
effects (MeshFeedback); the edge graph is fully connected, so HankinSolids no
longer needs it. `Segue::Dissolve` is deleted.

**Device result** (both build configs, two full 29-leg tours, 6,848 frames):
**0 spilled frames, 16 fps on all 18 holds and all 29 legs**. The morph is the
effect's lightest phase (legs 13.8–21.8 ms vs the 51.1 ms dissolve crossfade
it replaces); the worst frame is an ordinary hankin hold (59.0 ms). The
per-frame Conway op + `compile` costs **173 µs** — 0.3% of a frame, so §6's
headroom claim held and §4.1's "amortize construction across the preceding
hankin cycle" contingency is moot.

**Where the implementation overruled this design** (this document otherwise
describes the landed state):

1. **ADOPT on the two ambo-chain edges is non-destructive** (§2.1): the leg
   derives `ambo(seed)` at construction and `seed_base_` stays Platonic. A
   destructive ADOPT strands the walk on the return leg — no operator recovers
   cube/octahedron from a held cuboctahedron.
2. **Palette provenance is geometric, not emission-order** (§2.5): legs map
   provenance by departed-face centroids with a checked bijection. Emission
   order diverges from the held mesh's face order after dual-swap wandering,
   which misroutes colors and reads on screen as a *position jump*. Leg color
   targets key on the **arrival bookend** classification (hankin star-face
   classes), because the hankin classification can be coarser than the swept
   one (rhombicuboctahedron merges its 6+12 squares) and a per-slot forward
   mapping otherwise collapses both classes onto face 0's palette.
3. **`HS_PROFILE_ORDERED_CYCLE` walks a fixed 30-leg `ORDERED_TOUR`**, not
   `leg_index % m` arithmetic, which provably never reached node 15
   (rhombicosidodecahedron). The tour is static_asserted to cover all 18
   nodes, all five settle edges and every family bridge, and to wrap.
4. **`T_EPS_AMBO = 0.005`** clamps the truncate 0.5 end (no face degenerates
   there, and 0.02 leaves a multi-pixel gap); `T_EPS = 0.02` stays at the 0
   ends. Settle frames = 12 confirmed by measurement. `SNUB_BRIDGE_TWIST =
   −0.40` (settle rotation 23.4° → 17.3°); `SNUB_DODECAHEDRON_TWIST` stays 0
   to preserve §7.1 bitwise registry exactness (−0.25 would halve its residual
   rotation — an open cosmetic call, not a defect).
5. **Collapsed rosette faces are culled at draw time** (`SDF::Face`,
   |signed area| < 1e-5·radius²). §2.4's claim that the bookend swaps change
   no pixels was false as written: zero-area straps still rasterized as ~1 px
   lines along every edge. §7.6's bookend pin tolerates ±2 LSB per channel —
   the star 2n-gon's edge planes are fp-distinct from the base n-gon's, so
   bitwise pixel identity is unreachable.
6. **Settling legs run one eased clock** across sweep+settle (§2.3's separate
   linear settle left a visible near-stop at the seam).

## 1. Facts this design builds on

1. **Morphs happen between pure solids.** The interlace-angle sweep ends at
   angle 0, where every Hankin star point collapses onto its corner via the
   `p_corner` fallback (`core/mesh/hankin.h`, `update_hankin`): the hankin
   mesh at angle 0 *is* the base solid, rendered with midpoint vertices. The
   morph window therefore never needs the hankin layer — no per-frame
   `compile_hankin`, and the drawn mesh is a plain solid (≤ 92 faces —
   snubDodecahedron is the largest), far lighter than the interlace meshes
   already holding 16 fps.

2. **The registry is Conway chains, not baked meshes.** Every simple solid is
   `SolidBuilder(seed).op(params)` (`core/mesh/solids.h`). A morph that sweeps
   the final operator's parameter arrives at geometry produced by the *same
   code path* as the registry — for non-relax chains, bit-identical.

3. **Four operators are parameterized with geometric identity at t = 0:**
   `truncate(t)` (with `ambo == truncate(0.5)`), `expand(t)`, `chamfer(t)`,
   `snub(t, twist)`. Geometric identity only — at t = 0 they still emit the
   full expanded topology (duplicated vertices, zero-area orbit faces) with
   coincident positions, which is why sweeps start at T_EPS, never 0. No new
   operators are needed for v1.

4. **Topology is constant on the open interval.** For t ∈ (0, 0.5) a sweep has
   fixed connectivity; only vertex positions move. Face emission order is
   deterministic (primary faces first, in source-face order, then orbit/edge
   faces), and `MeshOps::compile`'s degenerate-face drop and orphan compaction
   are connectivity-driven, so face order and counts are identical every frame
   of a leg. Classification and palette assignment hoist to once per leg.

5. **`relax` preserves vertex count and order** (copy-then-relax-in-place), so
   a relaxed endpoint corresponds 1:1 to its unrelaxed form — the "settle"
   phase (§3.3) is a per-vertex slerp.

6. **Symmetry families connect only through tetrahedral shapes.** Conway
   operators preserve symmetry; the octahedral (order 24) and icosahedral
   (order 60) rotation groups have no containment relation, but chiral
   tetrahedral T (order 12) is a subgroup of both. The bridges:
   - `ambo(tetrahedron)` is exactly the regular octahedron (edge midpoints,
     normalized).
   - `snub(tetrahedron)` has icosahedron topology (20 triangles, 12 vertices);
     `relax` converges it to the *regular* icosahedron, same as the
     snubCube/snubDodecahedron chains already rely on relax for canonical form.
   - The `snub(tetrahedron, t, twist)` family also contains both bridge
     endpoints directly: the **exact** regular icosahedron (no relax) and, at
     (0.5, −π/3), the octahedron — the jitterbug bridge (§3).

7. **Existing utilities suffice:** `MeshOps::compile` (PolyMesh → drawable
   MeshState with face_offsets), `classify_faces_by_topology` PolyMesh
   overload, and the `Persist` arena-compaction idiom. No `compile_hankin`
   change is needed either: it emits exactly one star face per base face,
   first, in base-face order, which makes the bookend mapping the identity
   (§2.5).

## 2. Model

A **node** is a simple-registry solid, identified by canonical coordinates
`(seed, op, t)` — e.g. truncatedCube = (cube, truncate, 1/(2+√2)), cube =
(cube, identity, –). An **edge** is a single animated operator sweep between
two parameter values on one seed. The effect's persistent morph state is the
current node plus its **seed mesh** (a small PolyMesh, ≤ 120 vertices, in
`persistent_arena`) — kept in whatever orientation earlier legs produced, so
bridge arrivals need no realignment (the whole pipeline is
orientation-agnostic and the camera does a random walk anyway).

### 2.1 Reseed primitives (run in the leg-completion `.then()`)

- **ADOPT(node):** regenerate the clean endpoint mesh via its own chain, in
  the current seed's frame, into persistent arena as the new seed. Used when
  the next leg operates *on top of* the arrived solid (cuboctahedron →
  truncatedCuboctahedron) and at family bridges (arrived octahedron /
  icosahedron becomes the new family seed).
- **DUAL_SWAP:** at an ambo crossover, `seed := MeshOps::dual(seed)`,
  `t := 0.5`, continue sweeping down. `ambo(dual(seed)) == ambo(seed)`
  geometrically, so the crossover is seamless. Enables the classic dual
  traversal cube → cuboctahedron → octahedron as one continuous truncation.

### 2.2 Clean-endpoint swap

Sweeps are clamped to `[T_EPS, t_end]` (and away from exactly 0.5, where
`truncate` short-circuits to `ambo` with different emission order). At a leg
boundary the near-degenerate mesh (coincident vertex pairs, subpixel faces)
is swapped for the clean endpoint topology (seed mesh, or `ambo(seed)`, or
the registry-chain output). The swap is pixel-invisible by construction: the
degenerate elements it removes are below one pixel. `T_EPS = 0.02` initial;
tune against framebuffer dumps.

### 2.3 Settle phase (relax-terminated chains)

Five Archimedean chains end in `.relax(50)` (truncatedCuboctahedron,
snubCube, rhombicosidodecahedron, truncatedIcosidodecahedron,
snubDodecahedron), and both family bridges rely on relax for canonical form.
Relax is not sweepable, but by fact 5 it is slerpable: compute
`relaxed = relax(op(seed, t_end), 50)` **once** at leg construction, then over
the final `SETTLE_FRAMES` of the leg slerp each vertex from the unrelaxed
sweep position to its relaxed counterpart. Reverse legs un-settle first
(slerp relaxed → unrelaxed over the opening window), then sweep down.
`relax` never runs per frame.

Note snubDodecahedron's chirality comes entirely from relax (registry uses
`snub(0.5f)` with zero twist): with a plain settle the twist arrives late, as
a snap. Mitigation: sweep a nonzero cosmetic twist during the leg and let
settle close the (smaller) gap — a tuning choice, decided on framebuffer
dumps (§7.7).

### 2.4 Leg anatomy

```
[bookend-in]       at hankin-cycle end: angle forced to exactly 0, then
                   hankin mesh → node base mesh (= op(seed, t_node)), one
                   frame, pixel-identical; the hankin↔base palette mapping
                   (the identity, §2.5) applies here, owned by HankinSolids
[swap-in]          base mesh → swept solid at t = T_EPS; emission-order
                   identity on primary faces
[sweep, N frames]  per frame: op(seed, t(frame)) in scratch → compile →
                   draw with the leg's hoisted classification
[settle, S frames] per-vertex slerp to relaxed endpoint (if applicable)
[swap + reseed]    clean endpoint mesh; ADOPT / DUAL_SWAP as tabled
[bookend-out]      arrived base mesh, one frame; base→hankin mapping, then
                   compile_hankin + start_hankin_cycle(), first sample
                   forced to angle exactly 0
```

Initial values: N = 48, S = 12 (0 for non-relax legs). There is no dissolve
or crossfade: exactly one mesh is drawn per frame, and every boundary swap
exchanges geometrically matching meshes. (Swap-in shows T_EPS-sized corner
cuts appearing on the first frame — ~1–2 px at H = 144, reading as motion
onset; T_EPS is tuned against framebuffer dumps, §7.7.)

Two exactness/pacing notes. (1) Bookend frames **force** the interlace angle
to exactly 0 rather than trusting the sweep's boundary samples — the
sin_wave Mutation over 64 frames can land its final drawn frame at
~0.002 rad, off the flat `p_corner` branch, which would silently break swap
exactness; §7.6a pins this. (2) The cycle's rhythm changes: today's boundary
is a 16-frame MeshMorph, the leg is N + S = 60 frames against a 64-frame
hankin sweep — deliberate (the morph is now the show), and N is the knob if
it reads slow.

The bookend swaps change no pixels — the angle-0 hankin mesh and the base
mesh draw the same shapes — and the hankin mesh is display-only in this
design (nothing geometric ever consumes its collapsed star points or
zero-area fillers). What the bookends buy is factoring: ConwayMorph becomes
completely hankin-agnostic (it exchanges only base and swept meshes, mapped
by emission order), while the hankin↔base identity mapping lives entirely in
HankinSolids at cycle boundaries — each swap independently testable (§7.6).

### 2.5 Palette continuity (what the old dissolve was masking)

A boundary swap exchanges meshes whose positive-area faces are geometrically
identical, so the only discontinuity it can produce is color — a face
landing in a different class and therefore a different palette. Continuity
is exact if palettes follow face provenance across each swap:

- **Births are free.** New face families always enter at zero area — the
  vertex/edge faces of a Conway sweep at t = T_EPS, and the hankin filler
  faces at angle 0 — and grow from nothing. They take fresh shuffled palette
  indices with no pop. (Color turnover for *surviving* faces comes from the
  palette crossfade, §2.6.)
- **Bookend swaps (hankin ↔ base), at cycle boundaries.** The hankin faces
  with positive area at angle 0 fill their base faces 1:1 — and the mapping
  is the **identity**: `compile_hankin` emits exactly one star face per base
  face, first, in base-face order (faces 0..F−1), and the remaining rosette
  faces are zero-area births at angle 0. No `compile_hankin` change needed;
  the emission-order invariant is pinned by a unit test instead (§7.6a),
  with a per-face provenance array as the defensive fallback if that order
  ever has to change. Owned by HankinSolids; ConwayMorph never sees it.
- **Leg swaps (base ↔ swept, swept → swept).** Primary faces correspond by
  emission order (swap-in, swap-out, ADOPT) or by class signature at the
  shared solid (DUAL_SWAP at the ambo point, where side count + class size
  are unambiguous). Owned by ConwayMorph.

A class merge across a swap (two source classes feeding one destination
class) is not voted away: the mapping is kept **per face** (a from-palette
index per face, F ≤ 92 during a morph), so merged faces start from their own
colors and converge during the crossfade (§2.6). A split is exact (both
halves inherit the same palette).

### 2.6 Palette crossfade (color turnover inside the leg)

The provenance mapping (§2.5) anchors the *from*-state; the crossfade frees
the *to*-state. Each leg shuffles a fresh target palette assignment for its
classification, and over the sweep every face's color glides from its
inherited (mapped) palette to the target:

- **Blend weight** w(frame): 0 over the first ~20% of the sweep, smoothstep
  to 1 by ~80%, exactly 0 at swap-in and exactly 1 well before swap-out — so
  both boundary swaps remain exact and the color change hides mid-motion,
  where every edge on screen is already moving.
- **Mechanics**: not per-fragment. Once per frame, pre-blend the palette
  ramps for each distinct (from, to) pair into scratch (a handful of classes
  × a small ramp — hundreds of lerps, blended on the baked entries); faces
  shade from the blended ramp via the same single lookup as today. The
  fragment path is unchanged.
- **Per-face from-indices** absorb class merges (§2.5): merged faces simply
  reference different from-ramps until w reaches 1.
- **Newborn faces** open in the color of the face they are born inside: each
  newborn class inherits the from-palette of the departed face nearest its
  first member, so a T_EPS-wide birth is a sliver of the underlying color
  rather than a target-colored one.
- **Collapsing faces** are the mirror: a face with no counterpart at the
  closing bookend takes the *target* class of the arrival face it collapses
  onto (nearest arrival centroid among the surviving prefix), so it closes
  into its host's color. A T_EPS sliver is not zero area — on
  dodecahedron ← rhombicosidodecahedron the collapsing faces still cover
  ~12% of the lit pixels on the final frame — so a target class of their own
  would read as a hard color lattice that pops away at the swap.
- By swap-out the assignment has already landed on the leg's target, so the
  forward mapping into the new hankin cycle carries it verbatim.

## 3. Edge table (v1)

All edges bidirectional. Every intermediate topology equals an endpoint
topology already proven at 16 fps, so mid-sweep raster cost is bounded by
shapes in the existing 18/18 clean-hold set. The 23 edges cover all 18
simple-registry solids — every node is a registry entry and vice versa,
the completeness property the §7.8 soak relies on. The two `.bevel(t)`
registry chains decompose exactly (bevel = truncate ∘ ambo, conway.h), so
the ADOPT-then-truncate legs below are the registry code path, decomposed.

### Octahedral family (seed: cube, octahedron)

| Edge | Op sweep | Notes |
|---|---|---|
| cube ↔ truncatedCube | truncate 0 → 1/(2+√2) | |
| cube ↔ cuboctahedron | truncate 0 → 0.5 | arrive on clean `ambo(cube)` |
| truncatedCube ↔ cuboctahedron | truncate 1/(2+√2) → 0.5 | partial sweep, no degenerate end |
| cube ↔ rhombicuboctahedron | expand 0 → 2−√2 | |
| cube ↔ snubCube | snub (0,0) → (T_SNUB_CUBE, 0.28) | + settle |
| cuboctahedron ↔ truncatedCuboctahedron | truncate 0 → 1/(2+√2) | ADOPT cuboctahedron as seed; + settle |
| octahedron ↔ truncatedOctahedron | truncate 0 → 1/3 | |
| truncatedOctahedron ↔ cuboctahedron | truncate 1/3 → 0.5 | partial sweep, octahedron seed |
| octahedron ↔ cuboctahedron | truncate 0 → 0.5 | DUAL_SWAP cube ↔ octahedron at crossover |

### Icosahedral family (seed: dodecahedron, icosahedron)

| Edge | Op sweep | Notes |
|---|---|---|
| dodecahedron ↔ truncatedDodecahedron | truncate 0 → 1/(2+φ) | |
| dodecahedron ↔ icosidodecahedron | truncate 0 → 0.5 | |
| truncatedDodecahedron ↔ icosidodecahedron | partial truncate sweep | |
| icosahedron ↔ truncatedIcosahedron | truncate 0 → 1/3 | |
| truncatedIcosahedron ↔ icosidodecahedron | truncate 1/3 → 0.5 | partial sweep, icosahedron seed |
| icosahedron ↔ icosidodecahedron | truncate 0 → 0.5 | DUAL_SWAP dodecahedron ↔ icosahedron |
| dodecahedron ↔ rhombicosidodecahedron | expand 0 → 2−√2 | + settle |
| dodecahedron ↔ snubDodecahedron | snub (0,0) → (0.5, twist_c) | + settle; twist_c cosmetic, see §2.3 |
| icosidodecahedron ↔ truncatedIcosidodecahedron | truncate 0 → 1/(2+φ) | ADOPT icosidodecahedron; + settle |

### Tetrahedral family and bridges

| Edge | Op sweep | Notes |
|---|---|---|
| tetrahedron ↔ truncatedTetrahedron | truncate 0 → 1/3 | |
| tetrahedron ↔ octahedron | truncate 0 → 0.5 | **BRIDGE** tetra → octa; ADOPT octahedron |
| truncatedTetrahedron ↔ octahedron | partial truncate sweep | |
| tetrahedron ↔ icosahedron | snub (0,0) → (0.5, twist_b) | **BRIDGE** tetra → icosa; + settle (relax canonicalizes); ADOPT icosahedron |
| icosahedron ↔ octahedron | snub (t_jb, w_jb) → (0.5, −π/3) | **BRIDGE** (jitterbug) — see below; ADOPT at both arrivals |

**Jitterbug bridge (landed 2026-07-16)** — the third family bridge and the
first *direct* icosahedral ↔ octahedral crossing. The `snub(tetrahedron, t,
twist)` family is icosahedron-topology (V12/F20/E30) on the whole sweep and
contains both endpoints:

- The **exact regular icosahedron** — all 30 edges equal, no relax — at
  t = 0.7299092432622736, twist = −0.3881395153701886 (double-refined; pinned
  as `T_JITTERBUG_ICOSA` / `TWIST_JITTERBUG_ICOSA`). It matches the node's
  canonical `relax(snub)` held form per emission slot to 1.28e-3 chord
  (**0 px**), so the icosa-end swap needs no settle.
- The **octahedron** at (0.5, −π/3): the 12 vertices merge pairwise onto the
  6 octahedron vertices (cover 3.7e-8) and exactly the 12 edge-orbit faces go
  zero-area (hidden by the `SDF::Face` zero-area cull — endpoint render vs
  the clean octahedron measured **0 px**). Legs clamp at `T_EPS_JITTERBUG`
  (collapsing edge = 0.02 chord, T_EPS-sized), then clean-swap to the held
  octahedron; `ambo(tetra)` equals the registry octahedron exactly, so the
  octa endpoint is the held-octahedron invariant.
- ADOPT stores the octahedron at the octa arrival; the reverse (icosa)
  arrival regenerates the canonical `relax(snub(tetra, 0.5, −0.40))` held
  form, keeping the REGEN_TETRA invariant.

The cuboctahedron ↔ icosahedron jitterbug variant (twist-only snub on the
expand(T) form) stays dead: its quad → coplanar-triangle-pair swap creases a
√2 diagonal that does not shrink under any clamp — **61.7% pixel pop
measured** — the same defect class as the dropped pyritohedron cube leg.

Graph-walk policy: random edge from the current node, no immediate backtrack
**except at degree-1 nodes**, where an out-and-back is the only legal move
(six remain after the mirror edges: snubCube, rhombicuboctahedron,
truncatedCuboctahedron, snubDodecahedron, rhombicosidodecahedron,
truncatedIcosidodecahedron); bridge edges weighted so the walk changes
family every ~4–6 legs. Optional richness on the same patterns:
octahedron ↔ rhombicuboctahedron and icosahedron ↔ rhombicosidodecahedron
via expand, and expand(tetrahedron) → cuboctahedron as a further family
bridge — aesthetic calls, not gaps (they arrive off the registry code path,
so they'd test within tolerance like bridges, not exactly like §7.1).

## 4. Runtime architecture

### 4.1 `Animation::ConwayMorph` (core/animation/opleg.h, beside MeshMorph)

Constructor `(seed PolyMesh, EdgeSpec, Arena&, draw callbacks, frames,
easing)`:
- clones the seed into the arena (same survival contract as MeshMorph);
- runs the op once at the **clamped** arrival parameter — t_end pulled
  inside [T_EPS, 0.5 − T_EPS], never exact 0.5, where the ambo
  short-circuit changes emission order and face count and would hoist a
  classification misaligned with every swept frame's compile output — and
  classifies that PolyMesh (`classify_faces_by_topology`) into the arena:
  classified near arrival, not at T_EPS, so the 1°-angle-bucketed grouping
  matches what the viewer sees at leg end;
- if the edge settles: runs `relax(50)` once, stores the relaxed vertex
  array, and classifies the **relaxed** mesh instead — arrival geometry is
  the relaxed form, and for snubDodecahedron the unrelaxed (zero-twist) and
  relaxed (chiral) forms can bucket differently;
- builds the per-face from-palette indices via the leg-swap mapping (§2.5)
  and shuffles the leg's target palette assignment (§2.6). ConwayMorph deals
  only in base and swept meshes — it has no dependency on `CompiledHankin`;
  the hankin↔base bookends are HankinSolids' concern.

`step()` per frame: `ScratchScope` both scratch arenas → run the single op at
`t(frame)` (never the full chain — ADOPT guarantees one-op legs) →
settle-slerp if inside the window → `MeshOps::compile` into scratch → attach
the hoisted topology vector → pre-blend the (from, to) palette ramps at
w(frame) into scratch → invoke the draw callback. New profile tags:
`hk_conway_op`, `hk_conway_compile` (the existing `hk_draw_mesh` covers the
scan; the palette blend is a few hundred lerps and needs no tag).

`EdgeSpec` is a row of a `static constexpr` table: `{from_node, to_node,
seed_solid, op kind, t_from, t_to, twist_from, twist_to, settle, reseed}` —
23 rows, ~580 bytes, `.rodata` in flash on device.

**Storage decision — the table is flash, derived geometry is not.** Relaxed
endpoints, classifications, and provenance maps are computed at leg
construction, never baked: baked vertex data is frame-specific and ADOPT
deliberately leaves seeds in walk-dependent orientations; baking creates a
second source of truth against the §7.1 same-code-path guarantee; and the
construction cost equals today's `load_shape` spike (which runs the same
`relax(50)` chains) — already inside the 16 fps hold. If the §6 device
profile flags the spike anyway, the remedy is amortizing construction across
the preceding hankin cycle's frames, not flash baking.

### 4.2 HankinSolids integration

- `start_morph_cycle` selects a graph edge instead of a random registry index
  and constructs a ConwayMorph. `dissolve_`, `frame_tick_`'s dissolve salt,
  the two `MorphDrawFn` members, and the carousel's back slot become dead in
  this effect and are deleted (one mesh on screen means one draw path).
- New persistent member: `seed_base_` (PolyMesh) joins the `.then()`
  `Persist` compaction set alongside compiled_hankin / front slot / palettes.
- The completion `.then()`: clean-endpoint swap → reseed primitive → build
  the arrived base mesh → `compile_hankin` **from that mesh** (never a
  registry regenerate — orientation adoption after bridges) → classify →
  map palettes forward per §2.5 (the crossfade has already landed on the
  leg's target assignment, §2.6; fresh shuffles only for newborn classes) →
  compact → `start_hankin_cycle()`.
- The per-frame count `HS_CHECK` in `start_hankin_cycle` is untouched — the
  morph phase manages its own scratch-backed storage and never re-binds
  persistent arena.

### 4.3 Arena story

Per-frame swept meshes, HalfEdgeMesh, op scratch, `compile` output, and the
scan's FaceScratchBuffer all live in the scratch arenas under LIFO scopes;
peak is max-of-phases, not sum. Largest sweep allocation: truncate on
icosidodecahedron seed (V_out = 2E = 120, I_out = 360; the snub legs peak
at 60 V / 300 I). The per-solid scratch
high-water CI gate is extended with per-edge morph-frame peaks (§7.2).
Persistent additions: seed mesh + per-leg classification + relaxed vertex
array — all small and reclaimed by the existing compaction.

## 5. Operator-layer changes

Deliberately minimal:
- `T_EPS` clamp lives in ConwayMorph, not the operators; operators already
  accept the full parameter ranges (fail-fast guards stay as they are).
- Never evaluate `truncate` at exactly 0.5 mid-leg (ambo short-circuit
  changes emission order); the crossover uses the clean-swap + DUAL_SWAP.
- No kis changes (no simple-solid edge uses kis; the apex-height parameter is
  future work for Catalan legs, §9).
- `chamfer` stays unused (no simple-registry endpoint is a chamfered form).

## 6. Performance budget and gates

Added per-morph-frame work: one Conway op (O(I), I ≤ 360, HalfEdgeMesh build
included), `MeshOps::compile`, optional settle slerp — vector math measured
in microseconds against a mesh scan that renders a ≤ 92-face plain solid,
where today's morph frames render two *hankin* meshes (hundreds of faces)
under dissolve masks and hold 16 fps. Expected headroom is positive, but the
claim is gated, not assumed:

1. Device profile (teensy-profile regime) of the heaviest legs — dodecahedron
   → truncatedIcosidodecahedron (settle) and both snub legs — strict
   per-phase criterion: 🟢 only if every frame of every phase holds 16 fps.
2. Scratch high-water gate extended with per-edge morph-frame peaks.
3. Persistent-arena zero-growth check across a full graph walk (every node
   visited) — the soak test in §7.8.

## 7. Testing (all wired into tests/CMakeLists.txt + run_tests.cpp)

1. **Endpoint exactness.** For every non-bridge edge:
   `op(seed, t_end)` [+ `relax(50)`] equals the registry generator's output
   exactly (same code path, same seed frame) — vertex arrays, face_counts,
   faces.
2. **Topology-constancy sweep** (doubles as the degenerate-trap stress test).
   Per edge, ~16 samples of t ∈ [T_EPS, t_end]: constant V/F/I, closed
   manifold where the op requires it, all faces ≥ 3 sides, unit-length
   vertices, no traps.
3. **Settle correspondence pin.** `relax` output vertex order is the identity
   over its input — guards the slerp against future relax changes.
4. **Bridge convergence.** `snub(T, 0.5, twist_b).relax(50)` yields 12
   vertices / 20 triangles with all edge lengths equal within tolerance
   (regular icosahedron); `ambo(tetrahedron)` matches the regular octahedron
   within tolerance.
5. **Clean-swap invisibility.** `truncate(seed, 0.5 − ε)` vertices pairwise
   merge onto `ambo(seed)` vertices within tolerance; `op(seed, T_EPS)`
   primary faces match the seed's faces within tolerance.
6. **Boundary color continuity** — split along the §2.4 bookend factoring:
   (a) per **node**: unit-pin the identity-mapping invariant
   (`compile_hankin` emits star faces first, one per base face, in
   base-face order) and that the bookend frame's angle is exactly 0 (the
   flat `p_corner` branch — a final sweep sample at ~0.002 rad would break
   swap exactness); then framebuffer-diff hankin-at-0 vs base mesh —
   positive-area pixels match in color, not just coverage; (b) per **edge**:
   base vs op(seed, T_EPS) and the reseed swaps, same diff (exercises
   emission-order / class-signature mapping only); (c) unit: every mapping
   is total and deterministic, and the crossfade is exact at its endpoints —
   w = 0 reproduces the mapped from-state, w = 1 the target assignment.
7. **Visual keyframes.** Framebuffer-dump harness renders each edge at
   t = {ε, ¼, ½, ¾, end, settle-mid} to PNG — decisions about T_EPS, settle
   length, and snub twist cosmetics are made from renders, not theory.
8. **Device soak.** Full graph walk on Teensy visiting every node: no traps,
   per-phase 16 fps hold, zero persistent-arena growth. (Watch the
   stale-flash pitfalls: verify cycler markers, check log mtime.)

## 8. Risks

- **Settle aesthetics** on snubDodecahedron (late-arriving chirality) — has a
  planned mitigation (cosmetic twist sweep) and a decision procedure (§7.7).
- **Classification drift across a sweep**: angle buckets move with t;
  classifying at t_end fixes the grouping to arrival geometry, and the class
  *assignment* is frozen per leg, so no per-frame flicker is possible. The
  residual risk is a grouping at t_end that looks off near T_EPS — cosmetic
  and transient (those faces are near-degenerate there).
- **Emission-order dependence**: the bookend identity mapping leans on
  `compile_hankin`'s star-first, base-face-order emission — today an
  incidental property, not a contract; the §7.6 unit test promotes it to
  one. A per-face provenance array remains the defensive fallback if that
  order ever has to change.
- **Crossfade blend space**: ramps are blended on baked entries; if lerping
  baked colors muddies mid-blend hues on some palette pairs, blend earlier
  in the bake pipeline instead — an implementation detail, decided on
  framebuffer dumps like the other cosmetics.

## 9. Out of scope / future

- **Catalan solids**: kis apex-height parameter + dual legs (kis(X) is
  combinatorially dual(truncate(X))); HankinSolids cycles simple solids only.
- **Pyritohedron operator — DROPPED 2026-07-16 (owner), do not re-open.** The
  headline value (a *direct* cube → dodecahedron family crossing needing no
  Catalan) does not exist: the cube leg is measured-dead (below). What
  survives is one dodecahedron ↔ rhombicDodecahedron edge, which buys a
  crossing the tetra bridges already provide, at the cost of a Catalan node —
  and that visual territory is already shipping: `islamic_registry`'s
  `dodecahedron_hk72_ambo_dual_hk20` renders hankin on a
  rhombicTriacontahedron (`dual(ambo(dodecahedron))`) today. Evidence
  retained so the analysis is not redone. Vertices (±1,±1,±1);
  (0,±(1+h),±(1−h²)) + cyclic; 12 pentagons of 2 corners + 3 extras.
  - **h = 1/φ is the regular dodecahedron** (30 edges equal to 0.00e+00
    spread; edge multiset matches the registry to 2.65e-08) but in a
    cube-derived frame — arrival is bridge-tolerance, not §7.1-bitwise.
  - **h = 1 is exactly the engine's rhombicDodecahedron** (`dual(ambo(cube))`,
    8.5e-08, same axis-aligned frame, framebuffer diff **0 px**). The
    degeneration is pure pairwise vertex merge (12 extras → 6, ridge length
    2(1−h²) → 0) — the clean-swap class.
  - **h → 0 is NOT the cube.** Each cube face is tiled by two half-square
    pentagons with an internal ridge down the midline (a 90° great-circle arc
    per face, 6 total, plus collinear T-vertices). The crease is full-length
    at *every* small h — it does not shrink, so no T_EPS clamp makes it
    subpixel. Render at h=0.02 vs the plain cube: **91.6% of pixels differ**
    (dmean 63/255) — the pentagons are half-size faces, so `fragment_edge_dist`
    rescales the whole shading field, not just the crease. Same defect class
    as the dead cuboctahedron ↔ icosahedron jitterbug variant (§3), larger.
    **Do not re-attempt via T_EPS.**
  - Viable edge: **dodecahedron ↔ rhombicDodecahedron** (h: 1/φ → 1) — an
    icosahedral↔octahedral crossing, but it needs rhombicDodecahedron as a
    19th (degree-1) node, i.e. the Catalan-roster decision below, plus a new
    seed primitive (extract the inscribed cube frame from the held
    dodecahedron: the (±1,±1,±1) vertex subset, 5 valid choices); the landed
    jitterbug bridge (§3) already provides a direct crossing without either.
- **Hankin on Catalan solids is per-solid, not a blanket no** (measured
  2026-07-16, all 13 × 33 angles: zero traps, all vertices unit; the
  zero-fallback count was taken against a hard far-star cliff and has not been
  re-measured against the current blended, `plane_cross_sq`-gated fallback).
  Star-point chord ratio (cube control 1.52, dodecahedron 1.77):
  **rhombicDodecahedron 1.65 and rhombicTriacontahedron 1.80 are clean** and
  render as coherent interlaces across the full sweep. The kis family
  (1.5–3.2, non-monotonic) and the deltoidal/disdyakis/pentagonal snub-duals
  (peak **5.08**, well past the `STAR_FAR_RATIO_SQ` chord of 2.0 at which the
  fallback fully replaces the intersection) near-resonate mid-sweep — they
  render as sliver hairlines and star-face overlap. Mixed vertex degree is
  harmless; kis/snub-dual corner geometry is the discriminator.
- **Islamic registry**: chains end in hankin/multi-op stacks and would need
  per-frame hankin recompiles — explicitly rejected here.
