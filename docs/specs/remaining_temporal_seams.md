# Remaining temporal seams

Status: PROPOSED; inventory audited on 2026-10-17 against source
537f6606f. No additional effects are enabled by this spec.

## Scope and evidence

The implemented [pullback temporal sampling contract](pullback_temporal_sampling.md)
covers compatible composed effects, including GridSpace. This document records
the remaining coverage gaps and the design needed to close them on the physical
display and simulator together.

An effect without column-time sampling is susceptible to a sweep-boundary jump
when its image changes between frames. The inventory below establishes that
missing mechanism from source; it is not a claim that every preset has a
visually confirmed seam. Static images and slow or low-contrast motion can hide
the gap. No new device capture or visual seam survey was performed for this
inventory.

Temporal discontinuities must also be distinguished from spatial projection
singularities, aliasing, intentional discrete preset changes, and a snapshot
containing samples from different times. The simulator's visual response cannot
repair a discontinuity in the emitted image. Effect selection in the simulator
remains an atomic hard swap with replacement history.

## Outstanding inventory

| Effect | Source evidence | Required work |
|---|---|---|
| MobiusGrid | Composed effect with temporal sampling disabled. Animated Mobius coefficients use an interpolation overload that holds the previous lens until progress reaches one. | Separate temporal lens sampling from authored preset interpolation; validate the intermediate maps before opt-in. |
| SphericalHarmonics | One prepared harmonic frame, orientation, morph amount and amplitude are used by the entire shader draw. | Retain harmonic endpoints and resources; prepare the harmonic shader at each column time. |
| HyperLattice | One frame context and prepared trace state supply camera flight, rotations, parameters and palette to each frame's rays. | Retain motion and shading endpoints and prepare tracing per column for every supported trace path. |
| Voronoi | Sites advance once, then a frame-wide KD-tree and block candidate sets shade the image. | Evaluate moving sites at column time and make nearest-site candidate selection conservative over time. |
| Raymarch | Object centers, camera and volume orientations, spin and palette state are fixed during each volume draw. | Add column-time volume transforms and shading, with motion-conservative bounds. |

The source entry points are [MobiusGrid](../../effects/MobiusGrid.h),
[SphericalHarmonics](../../effects/SphericalHarmonics.h),
[HyperLattice](../../effects/HyperLattice.h), [Voronoi](../../effects/Voronoi.h)
and [Raymarch](../../effects/Raymarch.h).

The other effects outside the composed temporal path are BZReactionDiffusion,
Comets, DisplacementField, DreamBalls, Dynamo, Fishbowl, Galaxies, GnomonicStars,
GSReactionDiffusion, HankinSolids, HopfFibration, IslamicStars, MeshFeedback,
MindSplatter, MobiusRings, PetalFlow, RingShower, RingSpin, ShapeShifter and
Thrusters. They use forward geometry, particles, lattice simulation or feedback
paths. They remain candidates for temporal tearing, but cannot be fixed by
enabling the composed flag. A visual survey should rank their actual severity
before choosing a renderer-specific implementation.

## Shared design

Reuse the existing absolute-column timing function and optimized temporal
rotation preparation. For even width W, each opposing pair samples fraction
u = (x modulo (W / 2)) / (W / 2) between the previous and current endpoints.
At 288 columns the half-turn has 144 emission times: the last pair of one sweep
and the first pair of the next differ by one 1/144 frame interval. Clipping and
worker partitioning must never renormalize this time to their local width.

Expose a small shared retained-endpoint and column-preparation facility to
manual inverse renderers. Keep effect-specific state sampling as typed policies
outside the scan loop. Reuse `Scan::Shader::draw_columns` where its sampling
semantics fit; extend the generic scan entry points when cached directions or
spatial supersampling are needed. All spatial subsamples belonging to one
physical column share its emission time. Avoid one-off effect copies of timing,
quaternion preparation or resource-history ownership.

Advance simulation, timelines and random state exactly once per frame. Column
preparation is pure evaluation of retained state; it must not step animation,
consume random values, mutate a shared palette, or rebuild whole-frame resources.
The first frame primes both endpoints during initialization. Runtime fidelity-off
rendering uses the current endpoint while continuing to maintain history.

Retained endpoints must own mutable palette, noise and other lookup data that
can be rebaked before shading completes. Prefer sharing immutable resources and
existing endpoint color interpolation. Do not solve temporal discontinuities by
crossfading final rendered images: that duplicates moving edges rather than
evaluating their intermediate geometry.

Treat continuous motion separately from discrete state changes. Existing preset
semantics stay intact. Explicit resets, topology changes and reseeding need a
defined history reset or a supported transition; do not interpolate unrelated
object identities or invalid modes. Continuous coordinate wrapping must preserve
the field, including changes in its period.

## Effect adapters

### MobiusGrid

Add a temporal interpolation customization distinct from ordinary parameter
interpolation. Keep the existing authored lens snap policy. For the current
circular animation, retain and sample the animation phase and its defining
coefficients so intermediate lenses follow the actual path. This avoids replacing
a circular trajectory with a chord between coefficient snapshots.

A generic lens adapter must handle equivalent coefficient representations and
check nondegeneracy throughout its chosen path; valid endpoints alone do not
prove a valid interpolation. Use the existing determinant threshold and lens
limits as constraints. Unsupported transitions should explicitly reset history
or use an authored fade. Merely turning on temporal sampling while retaining the
held lens interpolation is insufficient.

### SphericalHarmonics

Retain orientation, the active harmonic pair, morph progress, amplitude and
palette resources. Prepare rotation and continuous controls once per column.
At a morph rollover, preserve the continuous field represented by the outgoing
pair at full progress and the incoming pair at zero progress; do not interpolate
integer harmonic indices. Manual mode snaps reset the history.

### HyperLattice

Retain unwrapped flight coordinates or explicit wrap deltas, six rotation
phases, parameters, preset gain and palette state. Sample periodic coordinates
only with a proven field period; do not assume every pattern shares the cubic
period. Prepare trace transforms once per column while preserving specialized
3D/4D, shell and experimental trace paths.

Crossing lists and other mutable scratch must belong to the column evaluation
or be proven invariant. Keep current per-ray step bounds, unfinished-ray
diagnostics and conservative culling. Begin with the supported shipping patterns;
gate experimental patterns separately until their timing and memory fit.

### Voronoi

Retain stable site identities, positions and rotation axes. Sample the known
rotation path per column rather than moving a frozen image. Nearest-two queries
must use positions at that time. A tree or block candidate set built only at
the current endpoint can miss the winning site at an intermediate time.

First build an exact small-site reference implementation. Then use a swept
spatial bound or a candidate scheme with an exact fallback, validated against
that reference. Rebuilding the complete KD-tree per column is a baseline for
correctness, not an assumed shipping solution. Reseeding or changing site count
resets history unless a stable identity mapping is explicitly supplied.

### Raymarch

Retain camera, object transforms, spin phases, continuous shape controls and
shading resources. Prepare each object's transform for the column before its
rays are marched. Volume culling must bound the whole interpolated interval, or
be evaluated for that column; endpoint-only bounds can remove visible motion.
Preserve object compositing order, antialiasing and shading coordinates.
Placement-solid and object-count changes require explicit history handling.

### Forward and simulation renderers

For forward geometry, retain identified primitive endpoints and evaluate
transforms at destination-column time. Clip rasterization to that column using
conservative swept bounds, or derive an inverse evaluator for suitable shapes.
Do not render the entire scene once for every column as the default design.

For simulation and feedback effects, retain consecutive simulation fields and
camera state and sample them at emission time where the display mapping allows
it. Measure the extra field storage first. Interpolating scalar simulation fields
can be useful, but does not establish correct particle lifetimes, topology or
advection. Those require their own temporal policy and reference tests.

## Acceptance and rollout

1. Record visible seam severity for each candidate with fidelity on, a fixed
   build and resolution, and deterministic motion. Include slow, fast and
   rollover cases; confirm important cases on the physical display.
2. Test sampled state against an independent time evaluator, including the last
   column pair of one frame and the first of the next. Use smooth reference
   fields to isolate timing from spatial image changes. Pin opposing-arm equality,
   first-frame initialization, clips crossing both half-turn boundaries, and
   single-engine versus segmented output.
3. Test palette rebakes, motion/clock wraps, morph rollovers, resets, discrete
   controls and fidelity toggles. Mutation tests must detect a held endpoint,
   local clip-time calculation and aliased mutable history.
4. For accelerated renderers, compare candidates and culling against an unculled
   exact temporal oracle across intermediate times. Include moving Voronoi
   winners and volumes outside their endpoint bounds.
5. Profile each adapter against the same optimized renderer with temporal
   sampling disabled, on the same device and configuration. Report incremental
   overhead, live mean/peak, spills, initialization cost, code size, persistent
   arena and stack usage. The [GridSpace device comparison](../profiles/shipping/profile_gridspace_teensy_2026-10-17.md) measured 11.2%
   incremental overhead; the earlier under-5% target remains unachieved and is
   not evidence that these adapters will meet it.
6. Land the shared manual-shader facility with SphericalHarmonics first; develop
   the independent lens policy for MobiusGrid; then tackle HyperLattice. Treat
   Voronoi and Raymarch as separate acceleration projects. Rank the forward and
   simulation backlog from the visual survey. Enable each effect only after its
   continuity, output parity, resource and device timing checks pass.
