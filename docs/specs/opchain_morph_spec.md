# Op-chain morph implementation note

**Status: LANDED.** IslamicStars builds authored solid recipes one operator at a
time. This document summarizes the shipped contracts; the source of truth is
`effects/IslamicStars.h`, `core/animation/opleg.h`, `core/mesh/recipe.h`, and
`core/mesh/solids.h`.

## Recipe model

`Solids::Recipe` stores a seed and authored `OpStep` values. The authored enum
contains `TRUNCATE`, `EXPAND`, `SNUB`, `CHAMFER`, `HANKIN`, `RELAX`, `KIS`,
`DUAL`, `AMBO`, `BEVEL`, `GYRO`, `META`, `NEEDLE`, and `ZIP`.
`expand_to_primitives()` lowers composite steps only when scheduling a build;
the stored recipe remains recognizable as the chain the author wrote.

`Solids::is_morphable_step()` is the admission gate. Unsupported or
out-of-range steps leave the recipe on the whole-shape fallback instead of
partially animating an invalid chain.

## Leg kinds

`Animation::OpLeg` has five construction-time kinds:

| Kind | Shipped use |
|---|---|
| `CONWAY_SWEEP` | Parameterized truncate, expand, snub, chamfer, and ambo-equivalent legs |
| `HANKIN_SWEEP` | Repositions compiled Hankin star points across the contact-angle sweep |
| `RELAX_SLERP` | Interpolates between identical-connectivity meshes around a relax step |
| `MEDIAL_SLERP` | Reconciles identity connectivity to exact authored endpoint positions |
| `GATED_SWAP` | Partition fallback for `kis` and `dual` |

Exactly one mesh is drawn per leg frame. Arrival topology, topology classes,
palette handoff state, and endpoint data are hoisted into arena-backed leg
state. Per-frame work produces the active mesh, pre-blends the distinct palette
pairs, and invokes the draw callback.

The palette bank has six entries. `MAX_BLEND_PAIRS` is
`PALETTES * PALETTES`, so the ramp table covers the complete 36-pair space;
only pairs actually used by a leg allocate scratch ramps.

## Truncate edge cases

The small-arrival and far-side recipe classes are both supported.
`T_EPS_TRUNCATE_FRAC` gives a `truncate(0.01)` leg a smaller positive birth so
it visibly sweeps instead of collapsing to a still image.
`T_EPS_TRUNCATE_FAR_MAX` bounds arrivals above 0.5, and
`truncate_off_pinch()` keeps a far-side sample on the constant-topology
truncate branch when it lands exactly on the ambo short-circuit.

An authored `bevel(0.5)` lowers through the exact ambo case rather than
emitting a topology-breaking truncate sample.

## Smooth kis/needle

The shipped smooth path replaces the visible partition swap for the supported
macro cases. A trailing `dual,kis` pair uses the dt bridge and a standalone
`kis` uses the dtd bridge. Each constructs identity connectivity, follows a
medial path, then runs a `MEDIAL_SLERP` reconcile leg onto the exact authored
endpoint positions. `needle` remains authored as a composite and lowers to
`DUAL,KIS`; the scheduler recognizes and spans that pair as one bridge.

The caller provides a checked nearest-vertex correspondence between the bridge
and authored endpoint. Connectivity remains fixed during reconciliation, so no
per-frame matching is needed. Palette provenance is carried through the same
arrival handoff used by ordinary legs.

## Renderer constraints retained from the design investigation

- Coplanar face partitions are not pixel-identical to their parent under the
  face-local antialiasing and source-over compositor. Interior partition edges
  remain visible even when shading gain is flattened.
- Pure radial motion does not provide a useful face-silhouette sweep under the
  gnomonic face representation; it changes shading inputs more readily than the
  projected boundary.
- `kis` and `dual` therefore need either the smooth bridge above or the explicit
  gated-swap fallback. A shading-only trick is not an exact replacement.
- Palette continuity is geometric. Departed-face centroids and arrival topology
  classes drive the checked mapping; emission order is not a stable identity.

## Historical dead ends

The following approaches were measured or structurally rejected and should not
be re-proposed without new renderer evidence:

1. Flattening gain to make a coplanar partition pixel-identical.
2. Animating a `kis` apex only by radial distance.
3. Parsing lossy registry names to reconstruct recipes.
4. Folding every relax step into its preceding leg despite vertex-count
   mismatches.

The earlier `truncate(0.01)` and far-side `truncate50d` failures are not on this
list: both were subsequently implemented as described above.

## Validation contract

Native coverage checks recipe replay against the authored generators, lowered
composites against their authored forms, topology and manifold invariants across
legs, palette-pair bounds, arena high-water limits, and end-to-end effect smoke.
Device profiling remains the authority for cadence and arena/ITCM acceptance.

The daydream solid editor owns its operation metadata and application logic in
`solid_codegen.js`; `solids.html` owns the UI and saved-chain state.
Those files, rather than historical line-number citations, are the browser-tool
sources of truth.
