# ShaderBall — unifying Liquid2D and Flyby

**Status: LANDED.** `effects/ShaderBall.h` is the sole ShaderBall implementation;
the earlier fixed stereographic implementation has been removed. The typed
pipeline carries 13 presets. Section 0 is authoritative for the authored
vocabulary, presets, and choreography, and for nothing beyond them: the shipping
renderer composes the reusable `Pullback::Pipeline` and operator catalog with
ShaderBall-owned `TopologyKey` and `ProgramDescriptor` selection, specified in
[the inverse-sampling pipeline spec](inverse_sampling_pipeline_spec.md). Nothing
in this file defines that architecture. Sections 1–13 preserve the original merge
and migration record; they are historical where they disagree with Section 0 or
the code. References to ShadierBall below name the former prototype from which
the final typed implementation was promoted. Section 11 remains the executed
roster checklist.

## 0. North star: authored field pipeline, pullback renderer

The composition engine, standard carriers, stage contracts, instrumentation
events, and reusable operator formulas live in `core/render/pullback.h`.
ShaderBall owns animation, frame preparation, resource lifetime, serialized
enums, topology admission, and its closed program manifest. Narrow empty state
providers are the only bridge from those effect-owned records into core
policies.

The north star is an authored effect pipeline, not the textual order of a
per-pixel shader. It should read the way an artist constructs the result:

```text
source animation → sample → warp → project → surface lens → colorize
                 → outer camera → viewer
```

The outer camera is last. The stack's implementation may traverse spatial and
deformation stages backward to answer a pixel lookup, but that traversal does
not redefine their authored order.

The fixed, typed stack must express every current ShaderBall preset without
turning ShaderBall into one opaque function slot. It must also retain
ShadierBall's additional functions, projections, and lenses. Continuous
parameters morph inside a topology; discrete slot tags choose the topology.

### 0.1 Migration baselines: both predecessors were pullback samplers

`Scan::Shader::draw` does not receive generated geometry or a forward-projected
sample. For each output pixel it supplies a world/view-space unit direction and
expects the final `Color4`. Both predecessor implementations therefore asked,
"which source coordinate and material land at this view direction?" The final
typed ShaderBall retains that pullback model.

For shipped ShaderBall, let `O` be the outer camera, `P` the current member
named `cam_inner`, `L` the glitch lens, and `project` the stereographic map. Its
exact spatial lookup is:

```text
outer_local = O⁻¹(view)
direct      = project(P⁻¹(outer_local))
lensed      = project(P⁻¹(L(outer_local)))
z           = lerp(direct, lensed, lens_mix)
```

It then evaluates:

```text
r_sq        = |z|²
(w, disp)   = noise_warp(z, r_sq)
field       = coupled_direct_function(w)
value       = pole_normalize(field, r_sq)
color       = colorize(value, disp)
```

The direct and lensed branches rejoin in projected-coordinate space, so warp,
function, value shaping, and colorization run once. This is the shipped lens
blend; the historical dual-pattern-sample description in §5 is superseded.

ShadierBall performs the subset:

```text
outer_local = O⁻¹(view)
z           = project_lens_blend(outer_local)
field       = selected_function(z)
pole_weight = pole_attenuation(|z|², pole_fade)
value       = (field * pole_weight + 1) / 2
coverage    = pole_weight²
color       = generated_palette(value)
color.alpha *= coverage
```

Its existing wander is consequently the outer camera. It is not an arbitrary
pre-lens animation that may be moved elsewhere when the framework grows.
ShadierBall currently applies this radial weighting and alpha coverage to all
three projection slots, including its bounded folded sinusoidal map. That
is shipped behavior even though the helper is named `pole_attenuation` and the
policy is geometrically motivated only by stereographic projection.

### 0.2 Authored order and lookup order

The authored stack describes field construction:

1. **Source animation** advances the source function's phase and orientation.
2. **Sample** evaluates the source function in its native domain.
3. **Warp** deforms that sampled field and emits deformation metadata.
4. **Project** maps the deformed field onto the sphere in an animated projection
   frame.
5. **Surface lens** distorts the projected spherical field.
6. **Colorize** converts the material sample and metadata into color and alpha.
7. **Outer camera** places the finished colored sphere in the viewer's frame.

As field operators, the authored composition is:

```text
Scene = OuterCamera ∘ Colorize ∘ SurfaceLens ∘ Project ∘ Warp
        ∘ Sample ∘ SourceAnimation
```

The hot loop evaluates a pullback of that composition:

```text
view lookup
  → inverse outer camera
  → surface-lens lookup
  → inverse animated projection frame
  → sphere-to-source projection lookup
  → warp lookup
  → source-animation lookup
  → sample
  → colorize
```

Equivalently, for one view direction:

```text
surface_dir = OuterCamera.lookup(view)
projected   = SurfaceLens.lookup_then_project(surface_dir)
source      = Warp.lookup(projected)
material    = Sample(SourceAnimation.state, source)
color       = Colorize(material)
```

Warp and projection therefore appear before `sample` in code even though they
appear after it in authored order. A warp is represented by its pullback map:
given a destination coordinate, return the source coordinate to sample. A
projection likewise maps a sphere lookup direction back to a source-domain
coordinate. The renderer never needs to forward-splat a dense intermediate
field.

Colorization is pointwise over the material sample and its carried metadata.
It is evaluated after the lookups have found that material, even though the
authored colored sphere is subsequently placed by the outer camera.

### 0.3 Retire the overloaded term "inner camera"

The proposed first phase, **source animation**, is not the current ShaderBall
member `cam_inner`. Source animation owns pattern travel, the independent second
pattern phase, and ShadierBall's function angle. It changes the field before it
is sampled.

The current `cam_inner` instead animates the frame in which the source field is
projected onto the sphere. In authored order it belongs to **Project**, inside
the surface lens and outer camera. In lookup order its inverse runs inside each
lens branch immediately before `project()`.

The north-star responsibility names are therefore:

- `SourceAnimation`: pre-sample clocks and source-domain orientation;
- `ProjectionFrame`: fixed alignment, spin, and optional projection-local
  wander represented today by `cam_inner`; and
- `OuterCamera`: the final authored sphere placement represented by
  `cam_outer` and by ShadierBall's current camera.

This vocabulary preserves the status quo while preventing "inner" from meaning
both source animation and a sphere-space camera.

These are phase labels, not a mandate for three new engine classes. Existing
`Timeline`, `Animation`, `Orientation`, and transformer facilities own the
corresponding state wherever their contracts fit.

### 0.4 Frame phasing

Every frame has one mutable preparation pass followed by an immutable lookup
pass:

1. **Choreography:** advance dwell/blend events, produce live continuous
   parameters, and latch completed discrete slot transitions.
2. **Source clocks:** integrate phase accumulators from those live rates and
   wrap each accumulator in its native domain.
3. **Spatial motion:** advance continuously running walks and prepare the
   `ProjectionFrame` and `OuterCamera` transforms.
4. **Resources:** step palette generation and prepare frame-constant noise,
   trigonometric, and transition state.
5. **Snapshot:** copy slot tags and hot values into an immutable `FrameState`.
6. **Pixels:** execute pullback lookups and pointwise sampling/colorization
   without mutating effect state.

Rates belong to presets; phase positions belong to the running effect. A preset
blend lerps rates and amounts, never accumulator positions. Entering a preset
changes velocity without resetting phase. Random walks likewise continue while
their gain is zero, so restoring gain does not cold-start a trajectory.

Manual parameter takeover pauses future preset selection. An in-flight preset
transition finishes, and source clocks, spatial walks, warp timers, and palette
resources continue.

### 0.5 Clock ownership

Named clocks replace ShadierBall's current assumption that every function can
share one travel phase and one angle:

| Clock | Domain | Owner | Current use |
|---|---:|---|---|
| `source_primary` | 2π | Source animation | both effects' main travel phase |
| `source_secondary` | 2π | Source animation | ShaderBall's independent drift phase |
| `source_angle` | 2π | Source animation | ShadierBall twin-wave/grid/spiral rotation |
| `source_noise_time` | basis-defined period | Source animation | noise-contour evolution and native wrap crossfade |
| `surface_noise_time` | 1 turn | Surface noise | sphere-space surface-noise evolution |
| `warp_outer_phase` | 1 turn | Planar Warp 1 | non-legacy warp animation |
| `warp_inner_phase` | 1 turn | Planar Warp 2 | non-legacy warp animation |
| `projection_spin` | 2π | Project | ShaderBall Y spin in the projection frame |
| `breathe_phase` | 2π | Colorize | ShaderBall breathe coordinate |
| palette resource state | provider-defined | Colorize resource | palette sequence, provider hue, and fade/dwell progress |

Slots consume only the clocks they declare. An unused clock may continue to
advance so a later transition into its consumer remains phase-continuous.
Every animated noise consumer names its clock and native period explicitly;
there is no implicit reuse of a `2π` source clock as a noise time axis.
The surface-noise slot owns `surface_noise_time`. Its rate is stored in
turns/frame, its position is captured in `FrameState`, and a discrete transition
forks and hands it off under the rules below. Because the clock drives a closed
loop through the noise domain rather than a linear time axis, it is exactly
periodic and needs no end-of-period crossfade. It may alias another noise clock
only when both the clock identity and complete `NoiseResourceKey` are explicitly
identical.
`source_angle` and `projection_spin` are deliberately distinct: the first
rotates a field in its native domain; the second changes how that field is
framed on the sphere. `breathe_phase` is likewise independent of palette
generation: ShaderBall's breathe rate is preset-controlled, while both effects'
palette providers and fade timelines advance on their own schedules and obey
their own pause policy.

### 0.6 Typed phase contracts

The authored phases exchange material meaning, but they do not each require a
new engine type. Most boundaries remain ordinary local values. The minimum
named records are:

```text
FrameState       { slots, live params, source phases, prepared transforms,
                   const resource bindings }
TransitionRuntime { mutable from/to clocks, orientations, mix, blend_mode }
TransitionFrame   { FrameState from, FrameState to, mix, blend_mode }
ProjectedLookup  { coords, region_id, component_id, boundary_flags,
                   fade_edge_distance, value_weight, flags }
StereoWarpResult { coords, displacement }       // existing
PlanarWarpResult { coords, path_length, flags }
MaterialSample   { value, coverage, path_length }
```

The pullback lookup path is:

1. `OuterCamera`: view direction → outer-local direction.
2. `SurfaceLens`: outer-local direction → active lookup branches.
3. `ProjectionFrame`: each branch → projection-local direction.
4. `Projection` kernel: each direction → complete branch `ProjectedLookup`.
5. Projection-specific surface-lens join:
   `join(ProjectedLookup a, ProjectedLookup b, mix) → ProjectedLookup`.
6. Post-join projection policy: recompute only scalar facts proven to be
   functions of joined coordinates.
7. `WarpProgram`: Planar Warp 1 lookup, then Planar Warp 2 lookup → source
   coordinate plus deformation metadata.
8. `Coordinate conditioning`: prepare bounded function arguments.
9. `Function`: source coordinate plus phases from `FrameState` → signed value.
10. `SignalWeight` + `ValueTransfer` + `CoveragePolicy`: signed value plus
    projection metadata → value and coverage.
11. `Colorizer`: value, coverage, and deformation metadata → `Color4`.

`coords` are pattern-plane coordinates, not necessarily a single continuous
Euclidean chart. `region_id` identifies the computational sector or
polyhedral face used by the formula; `component_id` identifies a connected
planar source component; `boundary_flags` identifies relevant `GLUED`, `CUT`,
`PERIODIC`, and `SINGULAR` ownership conditions; and `flags` records
projection-specific parity state. `fade_edge_distance` is distance to the
nearest fade-eligible `CUT` or `SINGULAR` boundary, never merely the nearest
glued edge. Projection-specific join logic owns any additional per-edge data it
needs. Region
identity alone neither proves nor disproves join compatibility. `value_weight`
is a finite projection-policy output in `[0, 1]`, not universal geometry.
`fade_edge_distance` is finite and nonnegative. Projection validation rejects
results outside either contract; Coverage does not silently clamp an invalid
projection result. All legacy ShadierBall slots initially
derive it from the joined coordinate's radius to preserve shipped output; new
bounded projections default it to 1 unless their own specification says
otherwise. Coverage does not belong to Projection: the same stereographic map
feeds both effects, but only ShadierBall squares the weight into alpha.
ShaderBall stereographic parity requires `r_sq = |joined.coords|²`, with that
same post-join radius feeding both noise-warp attenuation and `value_weight`.
Branch weights must never be interpolated in place of this finalization step.

Topology is never reconstructed from joined `Complex` coordinates. The
projection-specific join validates compatibility, interpolates coordinates,
and merges or selects region, component, boundary, parity, and fade-distance
metadata under that projection's exact rules. A fade distance may interpolate
only when both branches refer to the same boundary; otherwise the join is
incompatible. Post-join recomputation is limited to coordinate-derived scalars
such as stereographic `r_sq` and its radial weight.

The hot loop keeps the original `ProjectedLookup` beside `PlanarWarpResult`
rather than copying projection metadata through the warp program. Projection
metadata remains available until the value stage; warp-local flags may add fold
or singularity facts but never rewrite the originating projection region,
component, or edge. Coverage reaches final alpha, while `path_length` reaches
deformation-aware colorizers. A stage reports its displacement before
coordinate-addition rounding. The program accumulates those stage lengths
rather than recovering them by subtracting large rounded coordinates. With
legacy stereographic noise as the sole stage, `path_length` is the existing
helper's directly computed displacement scalar; the liquid colorizer consumes
that scalar without recomputing
`length(final_coords - original_coords)`. `MaterialSample` is an effect-local
merged shape, not a general engine type.

There are two clock models, selected by transition topology:

- A stable-topology parameter morph has one live state. Continuous rates lerp,
  and one authoritative clock integrates the resulting live rate, matching
  ShaderBall today. Clock positions never lerp.
- A discrete transition forks both endpoint clock positions from the
  authoritative state at transition start. Each branch advances at its own
  endpoint rate for the whole window, whichever endpoint the through-clear
  phase is currently drawing. At exact mix 1 the destination branch's
  positions become authoritative.

Spatial state follows the same distinction. A stable-topology morph integrates
one projection and outer orientation using the live lerped gains. A discrete
transition forks `inner_wander` and `outer_wander` orientations. The
underlying shadow `RandomWalk` animations still step exactly once and publish
one inner and one outer delta per frame; each endpoint integrates those shared
deltas with its own fixed endpoint gain. At exact mix 1 the destination
orientations become authoritative. A transform owner that cannot provide two
independently advanced endpoint states is not transition-admissible and the
discrete change snaps.

For the second model, mutable effect-local `TransitionRuntime` owns the forked
endpoint clock positions and integrated orientations that advance between
frames. Frame preparation derives an immutable `FrameState` for the endpoint
the current phase draws: slot tags, continuous params, prepared clock
values/transforms, and const resource bindings. Mutable resource
owners live outside both runtime branches and snapshots. Preset constants do
not drift underneath the transition. Each distinct mutable resource owner steps
exactly once per frame;
aliased endpoint bindings must not double-step the same `PaletteCycler`, noise
state, or animation. `mix` has exact 0 and 1 endpoints. The destination becomes
the sole active state only at exact 1, after which the source can be released.
Manual takeover edits the destination/live state and pauses future preset
selection; it does not partially rewrite the captured source endpoint.

### 0.7 Slots, parameters, and transitions

The north-star discrete slots are:

- `Function`: twin wave, rings, spiral, grid, and ShaderBall's coupled/direct
  field, plus noise contour and primitive lattice;
- `WarpProgram`: two fixed planar stages selected from the catalog in §0.7.1;
- `Projection`: folded sinusoidal, stereographic, gnomonic, Bonne, Peirce
  quincuncial, Airocean (Dymaxion), and equirectangular;
- `SurfaceLens`: none, glitch, twist, kaleidoscope, Möbius, and tangent noise;
- `SignalWeight`: none or projection weight;
- `ValueTransfer`: linear, ridge, iso-contour, or smooth bands;
- `CoveragePolicy`: opaque, projection weight, projection-weight squared, value
  cutout, or edge fade; and
- `Colorizer`: generated triadic, complementary, or analogous palette; cup,
  bell, ascending, or descending brightness envelope; and hue-shift mode.

Value policy owns distinctions projection cannot express. Given raw field
`f` in `[-1, 1]` and projection-supplied weight `w`, first compute
`s = clamp((f * signal_weight + 1) / 2, 0, 1)`, where `signal_weight` is `w`
for `PROJECTION` and 1 for `NONE`. Then apply one transfer:

```text
LINEAR:        value = s
RIDGE:         value = 1 - abs(2*s - 1)
ISO_CONTOUR:   value = 1 - smoothstep(width, 2*width, abs(s - level))
SMOOTH_BANDS:  value = 0.5 - 0.5*cos(2*pi*band_count*s + band_phase)
```

`width` has a positive floor, `band_count` is a positive integer,
`band_phase` is in radians, and all outputs clamp to `[0, 1]`. Coverage is
then:

```text
OPAQUE:                     coverage = 1
PROJECTION_WEIGHT:          coverage = w
PROJECTION_WEIGHT_SQUARED:  coverage = w²
VALUE_CUTOUT:               coverage = smoothstep(threshold-softness,
                                                  threshold+softness, value)
EDGE_FADE:                  coverage = smoothstep(0, edge_width,
                                                  fade_edge_distance)
```

`softness` has a positive floor. An `edge_width` of zero is a hard projected
edge; positive values use the smooth ramp above. ShaderBall maps exactly to
`PROJECTION + LINEAR + OPAQUE`; current ShadierBall maps to
`PROJECTION + LINEAR + PROJECTION_WEIGHT_SQUARED` for every legacy projection.
The former `UNWEIGHTED` behavior is `NONE + LINEAR + OPAQUE`. New bounded
projections normally supply `w = 1` unless an explicit edge treatment says
otherwise.

Colorize maps the transferred `value` to a palette coordinate independently of
coverage. The default `ASCENDING` envelope at frequency 1 is exact identity.
Other cases repeat the selected profile over the unit interval. For
`phase = wrap(min(value, 1 - ULP) * frequency)`:

```text
CUP:         palette_value = abs(2*phase - 1)
BELL:        palette_value = 1 - abs(2*phase - 1)
ASCENDING:   palette_value = phase
DESCENDING:  palette_value = 1 - phase
```

`frequency` is continuous in `[1, 32]`. Both the palette lookup and hue-rotation
lookup use `palette_value`; alpha coverage continues to use the material sample.

#### 0.7.1 Fixed projected-domain warp program

The projected-domain vocabulary is a fixed two-stage program, not an arbitrary
node graph. In authored order the stages are:

```text
source -> Planar Warp 2 -> Planar Warp 1 -> projection
```

The pullback shader evaluates the inverse order:

```text
ProjectedLookup.coords -> outer.lookup -> inner.lookup -> source function
```

`outer` is nearest Projection and `inner` is nearest Sample. Each position owns
one `WarpStageSpec { kind, basis, envelope, params, phase_rate, seed,
resource_id }`. Cost tier is derived from the selected kernels, bases,
integrators, wrap state, and transition shape; presets cannot self-declare a
cheaper tier. A zero amplitude is exact identity and bypasses all noise,
gradient, integration, and trigonometric work. Presets intending the fast path
store exact zero.

`resource_id` names an owner whose immutable `NoiseResourceKey` includes every
generator configuration field: basis/kernel type, seed, generator frequency
convention, channel offsets, octave/stencil version, and any mutable library
setting. Stage coordinate scale and clock value remain sample inputs only when
the kernel never pushes them into the generator object. Two stages or transition
endpoints may share an owner only when their complete keys match; otherwise they
require distinct simultaneously prepared owners. Alias compatibility, owner
count, and arena/flash cost are part of schema validation and transition
admission. No branch may call a configuration setter that clobbers another
branch's prepared resource.

The visual pipeline has compile-time `MAX_NOISE_RESOURCES = 9`, stored in a
fixed effect-owned array with no heap fallback. This covers the admitted maximum
of four distinct keyed owners per endpoint—Planar Warp 1, Planar Warp 2, noise
source, and surface noise—across a discrete transition, plus the hue-noise owner
whose key is endpoint-independent. Random-walk generators remain
separate dedicated owners and do not consume these slots. Before transition
capture, candidate validation deduplicates the union of endpoint keys and
rejects or snaps when it exceeds nine. All nine owner objects, bindings, and
any basis-specific state count toward RAM1/persistent-arena admission even when
the curated preset set normally uses fewer.

The effect-local results are:

```text
PlanarWarpStageResult { coords, delta, path_length, flags }
PlanarWarpResult      { coords, path_length, flags }
```

For original coordinate `p0`, evaluate `p1 = outer.lookup(p0)` and
`p2 = inner.lookup(p1)`. The final result has `coords = p2` and the sum of both
stage path lengths. Each `delta` is computed before adding it to that stage's
input, so it is not recovered by subtracting rounded coordinates. A direct
map's stage path is `length(delta)`; an integrated flow reports the sum of its
integration-segment lengths. Projection metadata remains the original
finalized post-lens `ProjectedLookup` beside this result.

The admitted stage catalog is deliberately orthogonal:

| Stage | Pullback definition | Continuous controls | Tier |
|---|---|---|---|
| `NONE` | identity | none | T0 |
| `AFFINE_FRAME` | `q = inverse(M) * (p - translation)` using a prepared nonsingular inverse | translation, rotation, log-scale x/y, bounded shear | T0 |
| `WAVE_SHEAR` | `q = p + A*sin(k*dot(d,p)+phase)*J*d` for unit `d` | amplitude, frequency, field angle, phase rate | T0 |
| `VORTEX` | rotate `p-center` by `2*pi*turns/(1+(r/radius)^2)` | center, radius, turns, center orbit/rate | T0/T1 |
| `VECTOR_NOISE` | add two decorrelated noise channels, then rotate the displacement vector | amplitude, scale, vector angle, time rate, envelope | T1/T2 |
| `CURL_FLOW` | integrate `dq/dtau = D*J*grad(N(q,t))`, `tau` in `[0,1]`, where `J(x,y)=(-y,x)` | signed flow distance `D`, scale, time rate, envelope | T2/T3 |
| `MIRROR_TILE` | rotate/offset, apply continuous triangular folds on x/y, rotate back | cell x/y, rotation, offset | T0 |
| `POLAR_CHART` | map to scaled radius and integral-harmonic angle; logarithmic mode uses `log(max(r, epsilon))` | radial scale, phases; positive integral harmonic and linear/log mode are discrete | T1 |

New controls use these closed ranges; writes outside them are rejected rather
than silently clamped into a different schema:

| Family | Normative ranges |
|---|---|
| affine | translation components `[-4,4]`, scale `[1/4,4]`, shear `[-3/4,3/4]`, rotation angle `[-pi/8,pi/8]` radians per phase cycle |
| wave shear | amplitude `[-4,4]`, angular frequency `[0,64]`, wrapped field angle/phase |
| vortex | center components `[-4,4]`, radius `[1/64,8]`, turns `[-4,4]` |
| vector noise | amplitude `[0,1]`, scale `[1/64,4]`, time rate `[-1/64,1/64]` turns/frame |
| curl flow | signed distance `[-1,1]`, scale `[1/64,2]`, time rate `[-1/64,1/64]` turns/frame plus the stability inequality below |
| mirror tile | cell dimensions `[1/64,8]`, offsets wrapped to their cells, wrapped rotation |
| polar chart | radial scale `[1/64,16]`, fixed radius floor `1/4096`, harmonic integer `[1,16]`, wrapped phases |
| noise contour | scale `[0,2]`, contrast `[0,8]`, time rate `[-1/64,1/64]` turns/frame |
| primitive lattice | cell scale `[1/64,8]`, shape blend `[0,1]`, softness `[1/1024,1]`, radius `[1/64,0.49]` in normalized cell units |
| value/coverage | levels and thresholds `[0,1]`, widths/softness `[1/1024,1/2]`, band count integer `[1,32]`, wrapped phase |

Möbius parameters use coefficient components no greater than 8 and must
satisfy a documented positive lower bound on `abs(a*d-b*c)` before that slot is
admitted. Authored presets use a canonical projective gauge: divide by positive
Frobenius norm, then rotate the common complex phase so the fixed-order
largest-magnitude coefficient is positive real; ties use coefficient order.
GUI edits retain their literal coefficients because a common scale cancels in
the projective transform; editing one coefficient never normalizes the others.
Coefficient lerp is forbidden. A one-lookup Möbius morph follows the prepared
matrix-group curve `M(t) = M0 * exp(t * log(inverse(M0) * M1))` under one fixed
logarithm branch and is admitted only when the complete curve preserves the
determinant and coordinate bounds. A branch change or failed proof uses the
through-clear transition or snap.

Every kernel input and intermediate must remain finite with each planar
component no larger than `WARP_COORD_LIMIT = 65536`. Every scaled FastNoiseLite
lattice coordinate, including central-difference offsets, must remain within
`NOISE_LATTICE_LIMIT = 2^20` before `floor`/integer conversion. Schema admission
proves these bounds from projection extent, stage order, and parameter ranges;
final source conditioning cannot repair an unsafe intermediate. No implicit
per-pixel saturation is allowed. A tuple whose bound cannot be proven is
invalid.

The curl sampler clamps each central-difference component to `[-G,G]`, where
`G = 4`. Curl with `n` integration intervals must satisfy
`scale * abs(D) * G / n <= 1/2`, and every predicted intermediate must also
remain inside `WARP_COORD_LIMIT`. This is an admission constraint in addition
to the slider ranges, not an adaptive runtime step reduction.

Transition admission covers the whole continuous path, not merely its two
valid endpoints. For every `t` in `[0,1]`, it proves the applicable determinant,
curl stability, warp-coordinate, and noise-lattice bounds. Curl uses the
conservative analytic bound
`max(scale0, scale1) * max(abs(D0), abs(D1)) * G / n <= 1/2` for a two-endpoint
morph, plus the path's coordinate bound. Operator-specific interval bounds
cover the remaining controls. If a continuous proof is unavailable, the
change is topology-incompatible and must use an admitted through-clear
transition or snap.

`AFFINE_FRAME` controls are authored as forward field placement, while frame
preparation stores the inverse pullback. Scale interpolates in log space with a
positive floor. Rotation is a signed angle in radians per phase cycle over
`[-2π, 2π]`, producing a continuous angular rate of `Speed × Rotation`; zero
leaves the frame fixed.
A singular or orientation-flipping affine endpoint is invalid.
`WAVE_SHEAR` fixes polarization perpendicular to field direction, making its
Jacobian determinant exactly 1 for every amplitude and frequency; a folding
wave displacement would be a different topology and is not implicit here.
Fold counts, polar harmonic, polar mode, noise basis, integrator, envelope,
`band_count`, and stage kind are discrete topology.

Noise stages select one fixed basis:

| Basis | Definition | Purpose |
|---|---|---|
| `SIMPLEX` | one OpenSimplex2 octave | broad smooth motion and the lowest-cost noise path |
| `FBM3` | three OpenSimplex2 octaves with fixed lacunarity 2, gain 1/2, and normalization by `1 + 1/2 + 1/4` | multiscale turbulence |
| `RIDGED3` | normalized three-octave `1-abs(noise)` spectrum with the same fixed octave weights | filament and ridge structure |

Channel offsets, derivative stencil, octave constants, and normalization are
kernel constants rather than sliders. Cellular is unavailable in the normal
firmware: `FASTNOISELITE_ONLY_OPENSIMPLEX2` makes runtime `SetNoiseType`
selection continue to sample Simplex. A future cellular basis is admitted only
by either removing that restriction build-wide with full ELF/cycle validation
or adding a standalone cellular kernel/type. Compiling incompatible
FastNoiseLite class definitions in different translation units is forbidden.
`VECTOR_NOISE` is a compressible displacement; `CURL_FLOW` is the
divergence-free field family and is not an alias for it.

Every basis publishes a canonical scalar range, but vector-noise channels also
have an exact zero-DC contract. Signed Simplex/FBM channels use their normalized
signed form. The nonnegative ridged basis forms each vector component
from a normalized difference of two decorrelated samples of the same basis,
so changing basis cannot add a constant translation. The resulting component
range is `[-1,1]`; fixed offsets and scale factors are part of the oracle.
Scalar sources use their own documented remap and need not share this vector
channel construction. Curl may consume the raw scalar basis because its
gradient removes a constant offset.

`CURL_FLOW` uses one signed flow-distance control over a fixed unit integration
interval; separate amplitude and integration length would be mathematically
redundant. Its discrete integrator tier is one Euler step, two midpoint steps,
or four midpoint steps. Step count and gradient stencil never degrade at
runtime. The chosen stencil spacing and basis time period are fixed in the
kernel's oracle. Every animated noise basis is mathematically time-periodic or
crossfades its native time-wrap seam; resetting a float accumulator alone is
not acceptable.

All new animated bases use the 3D sampling convention
`GetNoise(scale*x, scale*y, native_time)`; they do not switch silently between
2D and 3D overloads. Curl uses a central-difference gradient in scaled noise
coordinates with fixed `h = 1/64`:

```text
dNdx = (N(s*x+h, s*y, t) - N(s*x-h, s*y, t)) / (2*h)
dNdy = (N(s*x, s*y+h, t) - N(s*x, s*y-h, t)) / (2*h)
```

One gradient therefore costs four scalar-basis evaluations. One Euler interval
uses one gradient. Each midpoint interval uses two gradients, so the two- and
four-interval midpoint modes use four and eight gradients respectively.
`SIMPLEX` costs one base `GetNoise` call per scalar basis; `FBM3` and `RIDGED3`
cost three. `VECTOR_NOISE` costs two scalar bases for signed Simplex/FBM and
four scalar bases for paired-difference ridged channels. A time-wrap crossfade
multiplies the affected basis calls by two at its peak. Two stages add their
costs; a through-clear transition costs whichever endpoint it is drawing, so
its peak is the more expensive of the two.

Performance admission comes from cycle captures on the target device. Static
point estimates do not reject GUI combinations or rewrite selected stages.
Synthetic captures may still force time-wrap crossfades and transition pairs
when those paths need separate characterization.

Every non-legacy stage selects one envelope evaluated from the original
`ProjectedLookup`:

- `FLAT`: unit amplitude;
- `PROJECTION_WEIGHT`: multiply by the original `value_weight`; or
- `EDGE_FADE`: multiply by a clamped ramp from the original
  `fade_edge_distance`.

Progressive warp coordinates never move the projection pole or cut policy.
Warp-local folds and singularities append flags to `PlanarWarpResult`; they do
not replace projection region/component/edge metadata.

Stereographic `pole_fade` has one owner: the shared `StereoRadialPolicy` in
projection/value parameters. After the lens join it computes `value_weight`
from `r_sq`. No warp stage carries an independent pole-fade parameter; a stage
that attenuates by projection weight reads that same shared value through the
`PROJECTION_WEIGHT` envelope.

The initial source vocabulary adds two pure planar functions:

| Source | Definition | Controls |
|---|---|---|
| `NOISE_CONTOUR` | selected scalar noise basis remapped to `[-1,1]` | basis, scale, time rate, contrast |
| `PRIMITIVE_LATTICE` | repeated signed circle/box distance shaped to `[-1,1]` | circle/box blend, cell scale, radius, softness |

These source definitions are normative. For `NOISE_CONTOUR`, let
`n = clamp(N_basis(scale*x, scale*y, source_noise_time), -1, 1)`. Its output is
`n * (1 + contrast) / (1 + contrast*abs(n))`; this is exactly linear when
`contrast == 0`, remains odd, preserves `-1`, `0`, and `1`, and stays in the
declared range.

For `PRIMITIVE_LATTICE`, each axis uses the half-open repeated coordinate
`u = fract(cell_scale*p + 1/2) - 1/2`, so `u` lies in `[-1/2, 1/2)`. With one
shared half-extent/radius `r`, define:

```text
d_circle = length(u) - r
b = abs(u) - (r, r)
d_square = length(max(b, 0)) + min(max(b.x, b.y), 0)
d = lerp(d_circle, d_square, shape_blend)
value = 1 - 2*smoothstep(-softness, softness, d)
```

`shape_blend` is defined only on `[0,1]`; the square uses equal x/y half-
dimensions. Radius remains inside the half-cell while softness may span a cell
seam. The half-open ownership rule and the even SDFs give identical values on
neighboring cell seams; derivative parity may flip and is reported by the source
trait rather than claimed smooth.

Every source declares Cartesian-axis use, x/y periods or nonperiodicity,
rotation equivariance, polar-angle compatibility, and argument-conditioning
policy. `POLAR_CHART` is admitted only when the consuming source is periodic on
its angular axis and the integral harmonic lands its seam on an exact source
period.

The added sphere-space lenses reuse existing kernels where their direction is
valid: `MOBIUS` uses `mobius_transform` explicitly as a pullback. One persistent
lens does not justify a capacity-one pool. Ripple/burst lenses use
`TransformerPool` only when they actually spawn and reclaim independent animated
entities.

Sphere-space noise is a slot of its own rather than a lens, ordered against the
surface lens by `SurfaceNoisePlacement` (`BEFORE_LENS` or `AFTER_LENS`). Each
kind displaces the looked-up direction by an exponential map along a tangent
vector, so the result stays on the unit sphere:

| Slot | Pullback definition | Continuous controls |
|---|---|---|
| `NONE` | identity | none |
| `DIRECT` | step along a tangent sampled from the selected basis and rotated by `direction` | basis, scale, strength, rate, direction |
| `CURL` | step along the surface curl field, integrated by `EULER`, `MIDPOINT`, or `MIDPOINT_2X` | basis, scale, strength, rate, integrator |

Both kinds read `surface_noise_time`, which traverses a closed loop of radius
`NOISE_LOOP_RADIUS` through the noise domain. Scale is `[1/64, 8]`, rate is
`[-1/64, 1/64]` turns/frame, and `direction` is `[0, 1]` turns. Strength is
slot-dependent: `[0, 1/2]` for `DIRECT`, and `[-1/2, 1/2]` for `CURL`, whose
sign selects the flow direction. Zero strength is exactly identity and is the
only admitted way to hold the slot inert while its kind stays selected.

Composition rules are:

- field evaluation always uses the progressively warped coordinate;
- projection weight and edge fade always use the original `ProjectedLookup`;
- a same-stage continuous morph evaluates once with live parameters only after
  whole-path transition admission proves every intermediate valid; valid
  endpoints alone are insufficient;
- a discrete stage change may blend stage outputs and run unchanged downstream
  work once only when stage topology, fold/singularity schema, projection glue,
  source periods and metadata semantics are compatible;
- "unchanged downstream" means the complete downstream `FrameState` slice is
  branch-identical: slot tags, continuous params, clock values and rates,
  transforms, resource bindings, source, value/coverage policy, and colorizer.
  An operator-specific early join with unequal downstream state requires a
  separately proven commutation rule; otherwise both branches remain separate
  through their last differing consumer and normally blend rendered output;
- such an early join uses `coords = lerp(q_from, q_to, mix)` and
  `stage_path = lerp(path_from, path_to, mix)`; warp flags must be identical for
  that pixel because flags never interpolate. The unchanged stage then runs and
  total path is the sum of joined and unchanged stage paths. A flag mismatch
  uses the through-clear fallback or snap;
- otherwise the admitted through-clear fallback or exact snap applies; and
- if an outer-stage blend would cross an incompatible boundary consumed by the
  unchanged inner stage, the early join is not compatible.

Projection/warp/source compatibility is an explicit admitted table, not their
Cartesian product. Analytic stages require finite conditioned coordinates.
Noise across a `GLUED` edge is safe only when the field is equivariant under
that edge's coordinate-frame transform; the default is unproven. `CUT` and
disconnected components never coordinate-join. `MIRROR_TILE` reports derivative
parity flips. `POLAR_CHART` reports its origin singularity. Curl across a
reflected face requires explicit handedness correction because `J*grad` changes
parity. Conditioning runs after the complete warp program, and no NaN/Inf may
reach a source.

Each kernel ships with a double-precision oracle, exact identity test, finite
fuzzing over declared bounds, deterministic seed/phase/time-wrap tests,
composition-order and deformation-accounting tests, seam-neighbor tests,
compatibility-table tests, exact transition endpoints, and full device/ELF/
arena gates. Curl additionally proves gradient and integration output against a
double-precision reference. Diagnostic visual presets cover mirror-through-
vortex, vector-noise lattice, curl-flow grid, periodic polar source, and Möbius
over unwarped contours.

Out of scope are arbitrary graphs or equations, more than two planar stages,
bytecode or per-pixel type erasure, buffered feedback or previous-frame
advection, automatic Jacobians/lighting, silent octave or integrator reduction,
and implementing the entire catalog in one landing. Each kernel lands and is
admitted independently.

`ProjectionFrame` and `OuterCamera` are fixed positions, not responsibilities
of a function or lens. ShaderBall presets bind the same `wander` amount to
projection-local and outer walks for parity; separate gains may be exposed later
without changing phase order.

Topology controls never register GUI/deep-link writes directly against active
slot tags. They bind to an effect-local `RequestedConfig` shadow plus a dirty
set. At the next frame boundary the effect constructs one complete candidate
schema, validates all cross-field compatibility and resources, consults
transition admission, captures the active source state, and then either begins
the transition, commits an exact snap, or rejects the request without modifying
live state. Preset selection and deep-link import call the same whole-schema
apply path; they do not simulate an atomic change through sequential writes to
active fields. Continuous parameter writes may retain the existing direct
registered-pointer path only when they cannot invalidate the active schema.

Continuous parameters are grouped by their consuming phase. Presets store the
union required by compiled slots, while the GUI exposes parameters relevant to
the active topology. ShaderBall's `complexity` and `pattern_mix` remain inputs
to one coupled/direct function slot. Liquid and grid are not separate enums, so
current ShaderBall transitions stay within one topology.

Continuous angular parameters use canonical storage plus shortest-arc
interpolation with exact endpoint landing. An unwrapped accumulator may instead
be retained as the authoritative state and wrapped only in the immutable frame
value. Ordinary scalar `hs::lerp` across a periodic boundary is not valid.

Slot validity is explicit rather than assumed from the Cartesian product. The
compiled admission table covers projection, both warp stages, source traits,
value/coverage requirements, lens join, resource tier, and transition edge.
At minimum, `NONE` is valid for every projection. A stage whose definition
consumes stereographic radial attenuation is admitted only under stereographic
projection and must not be silently applied to Bonne, Peirce, Airocean, folded
sinusoidal, or gnomonic coordinates.
Preset tables are validated at compile time where possible; invalid GUI
combinations are rejected without changing the live state.

Enum values never lerp. A topology-changing transition blends at the earliest
boundary proven compatible:

- compatible lens branches may join in projected-coordinate space only when a
  projection-specific `join_compatible(a, b)` predicate proves that the
  straight coordinate interpolation does not cross a cut, singular boundary,
  incompatible orientation/parity frame, or disconnected component;
- functions over the same projected and warped coordinate contract may blend
  signed material values;
- compatible colorizers may use the premultiplied output blend defined below;
  and
- incompatible topologies fall back to a through-clear transition: the source
  fades to a fully cleared midpoint frame and the destination fades back up,
  so exactly one lookup renders per frame and the two endpoints are never on
  screen together.

Every transition and early join has direct exact endpoint branches before any
second lookup or blend math: `mix == 0` evaluates and returns only the source;
`mix == 1` evaluates and returns only the destination. This preserves bit-exact
endpoint color, avoids premultiply/unpremultiply round trips, and prevents an
unused expensive branch from running. Tests pin both output and kernel/noise
call counts at zero and one.

The fallback costs one shader evaluation per frame, not two, and its window
length is rounded up to an even frame count so the cleared midpoint lands on a
frame. Only one endpoint's persistent resources are prepared
at a time — the source's until the midpoint, the destination's after it.
Automated choreography may contain only admitted edges. An incompatible GUI
slot change whose pair is not admitted performs an exact single-frame snap to
the destination; it never substitutes coordinate interpolation. Lower
cadence/resolution or a cached-endpoint approximation requires its own future
specification and is not an implicit fallback.
Frequent or long transitions should prefer continuous amounts inside a stable
topology.

Legacy ShadierBall behavior is explicitly grandfathered: folded sinusoidal,
stereographic, and `SHADIERBALL_GNOMONIC_LEGACY` use
`LEGACY_PROJECTED_LERP`, which unconditionally interpolates the direct and
lensed planar coordinates exactly as shipped, ignoring fold and hemisphere
metadata. Strict `join_compatible` applies to new projections. Tightening a
legacy join later is a visual migration with new golden tests, not parity work.

The compatible-colorizer output blend is defined in premultiplied linear color,
even though the shader API returns straight-alpha `Color4`. Given endpoints
`(rgb0, a0)` and `(rgb1, a1)`, blend `rgb0*a0` with `rgb1*a1` and blend alpha
with the same mix; if the result must return through `Color4`, unpremultiply
when alpha is nonzero so `Scan::Shader::draw`'s final premultiplication
reconstructs the blended output. A zero-alpha result uses zero RGB. Straight
independent `Color4::lerp` is not the rule. Tests include unequal-alpha and
exact endpoint cases.

### 0.8 Projection family contract

Every projection's mathematical kernel implements the cartographic forward
lookup required by the pullback renderer:

```text
projection-local Vector -> Complex
```

This is the existing `stereo()`/`gnomonic()` convention: `Complex.re` and
`Complex.im` are planar x/y coordinates. Bonne and Airocean do not require
complex-number algebra, but `Complex` remains the engine's established compact
2D coordinate type. The projection slot wraps the coordinate only when it must
attach topology information:

```text
Vector -> Complex                         // pure projection kernel
Vector -> ProjectedLookup{coords, ...}    // slot result with metadata
```

"Forward" here means sphere to pattern plane. It does not contradict the
inverse rendering order: the shader has already pulled the output pixel back
through `OuterCamera`, `SurfaceLens`, and `ProjectionFrame` before making this
lookup. A projection need not provide plane-to-sphere inversion for rendering,
although a host reference may provide it for round-trip tests.

Each slot declares orthogonal topology traits rather than one lossy topology
enum: bounded/unbounded extent, periodic axes, fold/many-to-one behavior,
intentional cuts or interruptions, computational regions, connected source
components, and region orientation/parity. It must also define deterministic
ownership and classification of exact boundaries.

Projection kernels return their exact native unit-sphere coordinates: radius
or semimajor axis 1, native origin and offsets, x positive in the documented
eastward direction, and y in the documented northward direction. There is no
implicit normalization from an approximate or parameter-dependent full-map
extent. Existing projections retain their exact shipped coordinate scales; new
projection host oracles compare these native coordinates before any artistic
scaling. If comparable visual scale is desired, an explicit continuous
projection-coordinate scale is applied after the kernel and stored in the
preset. `pattern_freq` and coordinate conditioning follow that explicit scale
and any warp.

The projection roster is:

| Slot | Topology traits | Required controls | Boundary policy |
|---|---|---|---|
| `SINUSOIDAL` | bounded, folded, many-to-one | central meridian, legacy radial fade | preserve the shipped absolute-azimuth fold, pole collapse, and radial value/alpha policy |
| `EQUIRECTANGULAR` | bounded, one component with antimeridian cut | central meridian, legacy radial fade | expose the cut through metadata; longitude is unfolded and unscaled |
| `STEREOGRAPHIC` | unbounded, one component, singular pole | pole fade | preserve exact shipped mapping and radial attenuation |
| `GNOMONIC` | unbounded, antipodally folded, equator singularity | hemisphere policy, legacy radial fade | preserve the private ShadierBall kernel's exact divisor floor, component behavior, and radial value/alpha policy |
| `BONNE` | bounded, one component with antimeridian cut | signed standard parallel, central meridian | expose the cut through metadata; projection weight defaults to 1 |
| `PEIRCE_QUINCUNCIAL` | bounded, folded/tiled layout with oriented sectors | layout, central meridian, optional layout scroll | identify sectors/tiles, parity, and edges deterministically |
| `AIROCEAN` | bounded, interrupted polyhedral, multiple oriented regions | net orientation, central meridian | distinguish face/subface regions from connected components and interruption edges |

The shipped ShadierBall gnomonic kernel is not the core `gnomonic()` helper.
For parity it is provisionally named `SHADIERBALL_GNOMONIC_LEGACY`: it floors
the one shared signed `y` divisor to magnitude `1e-3`, divides both planar
components by that value, and performs no subsequent radial-sentinel clamp.
Replacing it with core `gnomonic()`, whose signed floor and radial sentinel
policy differ, is a separately reviewed visual migration.

Coordinate conditioning is part of the projection/source compatibility
contract, not an incidental helper call. All three legacy ShadierBall
projections multiply projected/warped components by `pattern_freq` and clamp
each function argument to `[-4096, +4096]`, preserving the shipped
`stereo_pattern_args()` behavior and fast-trig range bound. Each new projection
must declare its conditioning policy; it may reuse that bound, prove a tighter
native bound, or supply a different tested adapter.

#### Bonne

The Bonne slot is the spherical equal-area Bonne projection. Its
`standard_parallel` is signed and satisfies
`epsilon <= abs(standard_parallel) <= 90 degrees`; zero is invalid rather than
an implicit change to another projection. The reference preset value is
`+45 degrees`. `+/-90 degrees` are the corresponding Werner limits. Changing
the magnitude while staying on one side of zero is continuous. A sign change
or any transition through zero uses the incompatible-topology through-clear
transition.

The central meridian is a continuous angular control. The antimeridian cut is
intentional source topology, not numerical fallout; samples on it use a fixed
half-open ownership rule so neighboring pixels cannot flicker between sides.

#### Peirce quincuncial

The Peirce slot is the spherical conformal quincuncial projection. Its complete
world layouts are `DIAMOND`, `SQUARE`, `HORIZONTAL`, and `VERTICAL`; optional
hemisphere-only layouts may be added without creating new projection slots.
`SQUARE` is the reference layout because it exposes the quincuncial tiling
directly. Layout is discrete. Central meridian is continuous, while lateral
scroll is available only to layouts that define it.

The projection maps one hemisphere to the central square/diamond and divides
the other among the surrounding regions. Those folds are region boundaries.
The complete square or diamond tile is one connected planar component.
Reflected sectors receive distinct `region_id` and orientation/parity flags,
not fictitious disconnected chart IDs. Its projection-specific join predicate
may cross a `GLUED` or periodic sector edge when the coordinate frames agree;
crossing a cut, singularity, or incompatible reflection uses the
through-clear fallback.

Only the sphere-to-plane form is required in the device shader. This makes the
forward-only spherical reference implementation sufficient for rendering;
absence of a plane-to-sphere inverse is not a framework gap.

#### Dymaxion / Airocean

The canonical algorithmic identity is `AIROCEAN`; the artist-facing name may be
"Dymaxion / Airocean." This avoids treating every historical Dymaxion-like
icosahedral net as interchangeable. The slot follows the mathematically
specified Gray/PROJ Airocean map, including its fixed ocean-oriented cuts and
the subfaces used around Australia and Japan. A generic 20-face icosahedral
gnomonic map is not conformant.

`HORIZONTAL` and `VERTICAL` net orientations are discrete layouts. The shader
selects one face/subface, evaluates its local map, then places it in the fixed
net. Exact shared edges use a stable face-priority rule.
`fade_edge_distance` measures distance to the nearest fade-eligible
interruption or singular boundary, allowing value and coverage policy to
distinguish an intentional cut from a smooth interior. Planar source functions
still receive conditioned coordinates only; they do not inspect projection
metadata.

`region_id` records the selected face/subface, while `component_id` and the
per-edge `GLUED`/`CUT` classification describe the unfolded net. A face change
across a glued edge is not automatically incompatible; interpolation that
crosses a cut, singular boundary, or disconnected component is invalid because
it traverses unrelated source space. The projection-specific compatibility
predicate decides this from both branch results. All transitions between
Airocean and another projection use the general through-clear fallback.

Both sides of every paired cut use the authored edge-fade width. Bonne, Peirce,
and Airocean therefore meet with equal transparent margins rather than assigning
one side of the cut to pass beneath the other. Unpaired singular boundaries use
the same authored width.

#### Projection transitions, validation, and budget

Projection enum values and discrete layouts never interpolate coordinates.
Cross-projection transitions use the through-clear fallback from §0.7 when that
endpoint pair passes transition admission; otherwise the change snaps exactly.
Continuous
parameters may interpolate in one lookup while the topology and boundary schema
remain stable. Dynamic per-pixel region ownership may change as a central
meridian or other continuous parameter moves; the interpolated map resolves the
crossing with its deterministic boundary rule. A discrete layout, cut graph, or
projection change uses the through-clear fallback.

Each projection must ship with:

- host oracles against an authoritative double-precision spherical reference;
- named landmark, extent, central-meridian, and layout-orientation vectors;
- exact seam/face-edge tie-break tests plus finite-output fuzzing over the
  sphere;
- lens compatibility tests at region, glued-edge, and cut boundaries; and
- Teensy flash, ITCM, frame-cycle, and transition-cycle measurements.

The specification does not choose analytic code, approximation, or lookup
tables in advance. That choice follows error and device-budget measurement.
No new projection may silently degrade to another formula to meet budget.

Authoritative references:

- [PROJ: Bonne](https://proj.org/en/stable/operations/projections/bonne.html)
- [PROJ: Peirce Quincuncial](https://proj.org/en/stable/operations/projections/peirce_q.html)
- [PROJ: Airocean](https://proj.org/en/stable/operations/projections/airocean.html)
- [PROJ Airocean reference implementation](https://github.com/OSGeo/PROJ/blob/master/src/projections/airocean.cpp)

### 0.9 Retirement status

The typed implementation is now the sole ShaderBall. Retirement was accepted
on architectural and behavioral coverage rather than exact pixel parity with
the earlier fixed stereographic implementation. The final gate established:

- one authored/pullback phase model with outer camera last;
- named, continuously advancing source, warp, projection, breathe, walk, and
  palette clocks;
- the two-stage warp program, coupled/direct source, liquid colorizer, and all
  migrated looks in the preset bank;
- deterministic GUI edits and exact transition endpoints;
- premultiplied output blending for topology-changing preset transitions;
- projection topology metadata, host oracles, finite-output fuzzing, and seam
  ownership tests for Bonne, Peirce quincuncial, and Airocean;
- linear and squared projection-weight coverage, value cutout, and symmetric
  projection edge fade; and
- native effect, smoke, arena, stack, frame-budget, and WASM integration gates.

The earlier ShaderBall class and its dedicated white-box tests are removed.
The former ShadierBall class, registry name, test module, and simulator entry
are renamed to ShaderBall. Historical sections below retain the migration
record but do not define the shipping architecture.

### 0.10 Integration with existing engine architecture

The north star is assembled from existing engine mechanisms. A phase name is a
pipeline responsibility, not a request for a base class, virtual interface, or
manager object. Integration follows this ownership table:

| Responsibility | Existing owner | Integration decision |
|---|---|---|
| preset selection and continuous morphs | `Presets<>`, `Timeline`, `Animation::Lerp` | reuse unchanged; discrete slots use an admitted through-clear transition or exact snap |
| source, warp, spin, and color clocks | `Timeline` and existing `Animation` parameter drivers | reuse when the rate contract fits; keep explicit native-period integration for compound rates |
| projection and outer frames | `Orientation` and quaternion helpers | use `unorient()`/prepared conjugates directly in the pullback; do not add camera classes |
| animated sphere-space warps | `Transformer<>` and existing aliases such as `NoiseTransformer` and `MobiusWarpTransformer` | reuse when the effect needs spawned, composable animated entities |
| scalar fields on the sphere | `FieldTransformer<>` | reuse only for sources that are actually sums/compositions of spherical scalar fields |
| projection math | existing pure functions such as `stereo()` and the legacy private gnomonic kernel | extend the pure-function family; projection dispatch remains frame-constant |
| projected-plane warp | fixed effect-local `WarpProgram` over pure kernels such as `stereo_noise_warp()` | evaluate outer then inner; do not force `Complex` deformation through the `Vector -> Vector` transformer contract |
| source sampling | pure functions selected by the current function switch | retain; no generator object or per-pixel type erasure |
| palette motion and color lookup | `PaletteCycler`, palette recipes, and color utilities | retain separate persistent owners for liquid and generated-triadic state |
| procedural geometry construction | `hs::generate()` | not part of this shader pipeline |

#### Animation and frame preparation

`Timeline::step()` remains the sole timeline/event step. Existing
`Animation::Driver` is appropriate for a clock whose increment is one live
parameter times a constant. It is not a general expression driver: ShaderBall
clocks such as `speed * phase2_rate` and `speed * warp_time_scale`, and clocks
with periods other than `[0, 1)`, remain explicit accumulator updates in frame
preparation. They are wrapped in their native domains immediately afterward.
No `ClockBank`, `SourceAnimator`, or parallel scheduler is introduced merely to
move those few additions out of `draw_frame()`.

The required order aligns with current transformer lifecycle rules:

```text
Timeline::step(canvas)
  -> update/wrap explicit accumulators
  -> TransformerPool::prepare_frame() for pools whose params changed
  -> prepare Orientation conjugates and immutable frame locals
  -> Scan::Shader::draw(...)
```

The per-pixel lambda reads only prepared state. It never advances an animation,
calls `prepare_frame()`, mutates an `Orientation`, or allocates.

#### Resource ownership

Liquid and generated-triadic palette motion have two persistent
`PaletteCycler` owners outside `FrameState`, `TransitionRuntime`, and
`TransitionFrame`. This is
required because a cycler owns mutable provider sequence, fade/dwell counters,
generated palette slots, and display LUT state, and its pause binding is fixed
at initialization. The liquid owner binds `anims_paused`; the generated-triadic
owner does not. Every initialized owner receives exactly one step invocation in
the resource phase of every frame, even while no active endpoint references it;
its own pause policy may make that invocation a no-op. Aliased endpoint bindings
are deduplicated by owner identity, so no owner steps twice.

Frame snapshots carry only a resource identity and a const binding to the
owner's prepared palette/LUT. They never copy a `PaletteCycler`, mutable palette
storage, or gamut table. Non-owning `BakedPalette` copies do not create resource
independence and cannot stand in for two owners during a crossfade.

The persistent-memory gate includes at least both owners'
`generated_arena_bytes()` requirements, the ShaderBall gamut LUT, alignment,
and every other effect allocation that remains live concurrently. If those
resources do not fit, liquid/generated-triadic crossfade is not admitted and
the discrete change snaps under §0.7.

#### Where transformers fit

`Transformer<>` is specifically a fixed-capacity pool of independently spawned
animations whose hot path composes `Vector -> Vector` transforms in spawn
order. That is a direct fit for animated spherical surface distortion. An
effect may route a surface-lens slot to an existing transformer, for example a
Mobius or noise transformer, and the lookup renderer applies that transform
before projection.

Reuse does not erase pullback direction. `Transformer<>` itself only promises
function composition; it does not label a function as forward or inverse. When
bound into this shader, every supplied function is explicitly a destination-to-
source lookup map. Spawn order is consequently lookup composition order. For a
bijective authored transform, its inverse parameters/function must be supplied;
for a fold such as kaleidoscope, the lookup map is the effect definition because
there is no unique forward inverse. Existing transformer aliases may be reused
only when their mathematical direction matches that contract.

It is not automatically the right wrapper for every named Warp or Lens phase:

- fixed glitch, twist, and kaleidoscope lens slots remain direct pure
  `Vector -> Vector` functions unless they acquire spawned-entity semantics;
- `ProjectionFrame` and `OuterCamera` use `Orientation` directly because each
  is one persistent frame, not a pool of temporary warps; and
- stereographic noise warp remains a projected-plane operation returning both
  `Complex` coordinates and displacement. Adapting it through
  `Transformer<>` would change its domain, discard metadata, or add a
  sphere-plane round trip.

`FieldTransformer<>` has the same selective rule. It is appropriate when a
source is genuinely composed from multiple animated scalar fields evaluated on
the sphere. The current planar twin-wave, rings, spiral, grid, and
coupled/direct functions are one selected function each, so they remain inline
samplers rather than one-entity transformer pools.

#### Typed operator family

Architectural coherence comes from repeating the existing transformer pattern,
not from forcing every stage behind one interface. The engine already separates
an animation-owned parameter entity from the pure hot-path kernel that consumes
it. Preserve that split across domains:

| Operator kind | Mathematical shape | Shader role | Composition |
|---|---|---|---|
| spherical transformer | `Vector -> Vector` | surface lens / sphere-space lookup warp | existing `Transformer<>`, in lookup order |
| planar transformer | `Complex -> Complex` | one stage of the fixed projected-domain lookup program | direct pure kernel; optional sibling pool only if spawned entities are later required |
| projection kernel | `Vector -> Complex` | sphere-to-pattern-plane coordinate map | one selected pure map, not a transformer pool |
| projection slot | `Vector -> ProjectedLookup` | coordinate plus topology metadata | wraps one kernel; does not form a transform pool |
| planar source function | `Complex -> float` | pattern sampling | one selected pure function |
| spherical field | `Vector -> float` | scalar source/displacement field | existing `FieldTransformer<>` when multiple entities must reduce together |
| colorizer | material values -> `Color4` | terminal shading | one selected pure function/resource owner |

The common convention is:

```text
prepared Params + immutable input -> pure result
```

Animation changes `Params` before drawing; it does not live inside the sampling
call. Actual dispatch stays as frame-constant enum switches or compile-time
calls, but that does not require every arm to inline into the pixel closure. The
signatures above describe contracts, not stored function pointers or virtual
methods.

#### Hot-code placement and device admission

Only measured small kernels inline into the per-pixel closure. Ordinary
out-of-line code still occupies ITCM on Teensy, so every large projection,
noise, color, and transition leaf has an explicit placement admission:
measured shared ITCM, or `FLASHMEM`/equivalent XIP placement. Flash-resident
leaves report both cold-cache and steady-cache call-path cycles; they are not
assumed free merely because they save ITCM. An on-device A/B measurement may
promote a leaf only when total cycles improve and the full roster still fits.
Blanket `always_inline`/`HS_O3_FN` expansion of every slot arm is forbidden: a
frame-constant switch makes prediction cheap but does not remove inactive-arm
code from the instantiated closure.

The shipping bank retired former presets 5, 6, 7, 9, and 10. Two
source-matched selective-O3 captures run every
current preset green against the 62.5 ms display window, while the matched
global-O3 reference has five red presets and is rejected for shipping. The
shipping report preserves the former 17-preset indices and exact current
remap. The full-roster memory gate includes ShaderBall. Take the spill, peak,
and byte figures from the
[shipping profile](../profiles/shipping/profile_shaderball_teensy_2026-08-14.md),
[global-O3 reference](../profiles/O3/profile_shaderball_teensy_2026-08-14.md),
and [ITCM ledger](../ledgers/itcm_ledger.md) rather than restating them here.
The timing captures measure the authored bank directly; the ledger remains a
separate full-roster resource gate. Regenerate both after preset, kernel, or
roster changes and record the commit, build flags, roster, hardware, and capture
coverage beside the result.

Run a full-roster ELF/ITCM gate and on-device cycle profile after each
fixed-topology integration stage, each large kernel, and each new projection—not
only at the final retirement gate. A dedicated comparison build may alternate
the two implementations, but the merged replacement must pass the full shipping
roster ELF gate in ShaderBall's slot before removal. Transition admission uses
measured peak endpoint costs, not average frame time or curated-preset
intuition.

`TransformerPool<ParamsT, AnimT, CAPACITY>` remains the shared lifecycle base
for cases with independently spawned parameter entities. Its design already
supports sibling evaluation policies: `Transformer<>` folds vectors, while
`FieldTransformer<>` reduces scalar fields. If projected-domain warps later need
the same spawn/reclaim/composition behavior, the coherent extension is a
`PlaneTransformer<>` sibling derived from `TransformerPool`, not a generalized
render graph and not a conversion of plane coordinates back through the sphere.

That extension is deferred until independently spawned entity lifecycle is
required. The fixed two-stage `PlanarWarpProgram` is effect-local direct
dispatch and calls kernels such as `stereo_noise_warp()` itself. Two
preset-selected positions do not justify a transformer pool: creating two
capacity-one pools would add lifecycle and indirection without adding behavior.

A future `PlaneTransformer<>` must settle its reduction contract before it is
introduced: transform order, whether downstream metadata represents net
displacement or accumulated path length, and whether each transform evaluates
attenuation against the original projection lookup or the progressively warped
coordinate. Those choices affect rendered output and cannot be hidden in a
generic template default.

Projection deliberately does not derive from `TransformerPool`. Its kernel
changes domain from `Vector` to `Complex`, has no identity element in the target
type, and is selected once per lookup rather than composed as a sequence of
animated entities. The slot-level wrapper adds `ProjectedLookup` metadata only
where topology needs it. Projection follows the same pure-kernel/animated-
params convention while remaining a distinct typed map. Projection-specific
continuous parameters are ordinary preset values prepared by
`Timeline`/`Animation::Lerp`; changing the map or layout remains a discrete slot
transition.

Likewise, animated sampling means that `Timeline` advances source phases and
other sampler parameters before the frame snapshot. The sampler itself remains
a pure `Complex -> float` evaluation. Call it a source or field function in
code and documentation; calling it a generator would collide with
`hs::generate()`'s established procedural-geometry meaning.

#### Why `generate()` is not animated sampling

The engine's `hs::generate()` is a cold procedural-geometry wrapper. It resets
and scopes the shared scratch arenas and writes persistent generated output.
Calling it from a per-frame or per-pixel source sampler would violate both its
lifecycle and the shader's immutable hot-loop contract. It remains useful only
if a projection or colorizer needs a setup-time baked table; the resulting
table is then an immutable resource read by the shader.

FastNoiseLite is a noise generator in the ordinary sense, but it is not the
`generators.h` architecture. It remains effect- or params-owned state. Owners
must stay exclusive when configuration methods mutate the generator, matching
the existing RandomWalk and ShaderBall rule.

#### Abstraction threshold

The first integration keeps slot dispatch and orchestration effect-local and
reuses core math/helpers. It adds no virtual phase hierarchy, generic render
graph, stored `FunctionRef` chain, or generalized transformer domain. A shared
engine abstraction is justified only after at least two independent effects
need the same contract and the extracted form preserves hot-path inlining,
metadata, lifecycle ownership, and measured device timing.

The immediate architectural delta is therefore deliberately small:

1. one `ProjectedLookup` record for nontrivial projection topology;
2. pure projection functions and their frame-constant switch dispatch;
3. effect-local signal-weight, value-transfer, and coverage dispatch preserving
   the two shipped alpha behaviors;
4. effect-local mutable `TransitionRuntime` plus immutable
   `FrameState`/`TransitionFrame` snapshots, containing const bindings rather
   than copied resource owners; and
5. adapters to existing transformers only for actual sphere-space animated
   compositions.

## Historical ShaderBall merge design

Merge `Liquid2D` and `Flyby` into a single stereographic shader effect,
`ShaderBall`, whose presets span both effects' looks and the mixed looks
neither can reach today. The two effects already share their entire skeleton
(stereographic projection → noise warp → trig pattern → pole normalize →
palette); they diverge only along five axes, and every one of those axes can
be expressed as a **continuous float parameter**. That is the core design
decision: no mode enums, no per-preset function pointers — every preset is a
point in one parameter space, so the existing preset-lerp machinery morphs
between *any* two presets, including cross-family ones, with no discrete pop.

## 1. Why merge

- The classes are ~85% structurally identical (both note this in their own
  doc comments and share `stereo_noise_warp` / `stereo_pattern_args` /
  `pole_normalize_pattern` in core). Two copies of the skeleton cost flash,
  ITCM, test surface, and roster slots.
- The interesting space is *between* them: a glitch-lensed grid fly-through,
  a hue-warped cross-coupled liquid, a half-lensed slow orbit. Today those
  are unreachable; in ShaderBall they are just presets.
- One effect replaces two in every roster (HS_EFFECT_LIST, Phantasm playlist,
  GOLDEN_ROSTER, README gallery, daydream favorites), shrinking the
  maintenance checklist permanently.

## 2. Current anatomy

Shared skeleton (identical in both):

- `Effect(W, H, {.strobe = true})`, no pixel persistence; full-canvas
  `Scan::Shader::draw<W, H, 1>`.
- FastNoiseLite OpenSimplex2 warp source (but the generators' frequency,
  seeding, and sharing diverge — see the table below).
- Per pixel: project to stereographic plane → `stereo_noise_warp` →
  `stereo_pattern_args` + fast trig pattern → `pole_normalize_pattern` →
  palette lookup.
- Noise-time accumulator wrapped to `STEREO_NOISE_TIME_PERIOD`; trig phases
  wrapped to 2π (white-box tested in both).
- `Presets<Params, N>` bank + timeline-scheduled `Animation::Lerp`.
- Baked generative palette in the persistent arena.

Divergences:

| Axis | Liquid2D | Flyby |
|---|---|---|
| Camera | two `RandomWalk`s (inner + global), no spin | fixed base orientation + looping Y rotation (2π / 300 frames) |
| Lens | glitch lens (hemisphere fold + degree-3 rational map) between the walks | none |
| Warp noise generator | **shared with the inner walk**: `RandomWalk`'s ctor reseeds it (random seed) and sets its frequency; init then re-sets `SetFrequency(0.02)` | dedicated, untouched: FastNoiseLite defaults — frequency **0.01**, fixed default seed |
| Warp time | `t * 0.5` | `t * 0.3` |
| Pattern | cross-coupled: `sin(x + c·sin(y + φ₁)) · cos(y + c·cos(x − φ₂))` | separable grid: `sin(x + φ₁) · cos(y − φ₂)` |
| Phase rates | φ₁ += dt, φ₂ += 0.8·dt | φ₁ += speed, φ₂ += speed·drift (drift = live slider, default 0.7) |
| Color | `StaticPalette` + `BreatheModifier` (amp 0.15, cycle-phase driven) | raw palette + `alpha *= (1 − value)` + `hue_rotate(−displacement · hue_shift)` via gamut LUT |
| Palette recipe | STRAIGHT / COMPLEMENTARY / CUP / VIBRANT / 75 | STRAIGHT / SPLIT_COMPLEMENTARY / FLAT / MID / 42 |
| Preset choreography | repeating 90–150-frame timer; the 60-frame **staggered** lerp runs concurrently inside the period (actual hold between blends: 30–90 frames) | perpetual 480-frame **parallel** lerp loop, `.then()`-chained with zero gap |
| Persistent arena | palette only | palette + gamut LUT (131,072 B) |

## 3. The unifying observation: one pattern formula

Generalize the pattern to

```
pattern(p, φ₁, φ₂) = sin(p.re + c·sin(p.im + φ₁) + s·φ₁)
                   · cos(p.im + c·cos(p.re − φ₂) − s·φ₂)
```

with `c = complexity` (cross-coupling depth) and `s = phase_direct`
(direct phase feed). This is an **exact** superset:

- `c = complexity, s = 0` → Liquid2D's `sample()` verbatim.
- `c = 0, s = 1` → Flyby's `sample()` verbatim.

Both parameters are continuous and lerpable, so a preset transition morphs
the *function itself* smoothly. `stereo_pattern_args` still bounds `p`
first, exactly as today.

Phase advance unifies the same way: `φ₁ += speed`, `φ₂ += speed ·
phase2_rate` per frame (Liquid2D: `phase2_rate = 0.8`; Flyby: `phase2_rate
= drift = 0.7`). Warp time becomes `noise_time · warp_time_scale`
(0.5 vs 0.3).

## 4. Pipeline

Per frame (all frame-constant work hoisted out of the shader):

```
timeline.step()                            // walks, drivers, preset lerp
noise_time = fmodf(noise_time + speed, STEREO_NOISE_TIME_PERIOD)
sin_phase  = fmodf(sin_phase  + speed, TWO_PI_F)
phase2     = fmodf(phase2 + speed * phase2_rate, TWO_PI_F)
spin_phase = fmodf(spin_phase + spin_rate, TWO_PI_F)
cycle_phase = fmodf(cycle_phase, TWO_PI_F)  // Driver-fed, breathe consumer
inner = R_y(spin_phase) ∘ BASE ∘ slerp(I, inner_walk, wander)   // per frame
outer = slerp(I, outer_walk, wander)                            // per frame
```

Per pixel:

```
rv = outer.unorient(v)
sv = lens_blend(rv, lens_mix)          // §5
z  = stereo(inner.unorient(sv))
r_sq = |z|²
(w, disp) = stereo_noise_warp(z, r_sq, noise, warp_scale, warp_strength,
                              pole_fade, noise_time * warp_time_scale)
pattern = generalized formula (§3) on stereo_pattern_args(w, pattern_freq)
value = pole_normalize_pattern(pattern, r_sq, pole_fade)
value = min(value, 1.0f − ULP)         // §6: Wrap palette folds exact 1.0 → 0
c = palette_get(value)                 // §6: bank blend + breathe
c.alpha *= (1 − value * value_fade)
if (hue_shift ≠ 0 or disp-dependent) c = hue_rotate(c, −disp * hue_shift)
return c
```

Every per-pixel branch (`lens_mix == 0`, `complexity == 0`, `hue_shift ==
0`, `palette_pos` integral) tests a **frame-constant** value, so the M7
branch predictor makes the skips nearly free without duplicating the inner
loop into ITCM-costly specializations (see §9).

## 5. Camera and lens

Three continuous controls replace the discrete camera modality:

- **`spin_rate`** — Y-axis rotation rate in rad/frame. Flyby: 2π/300 ≈
  0.0209; Liquid2D: 0. The rotation is an explicit wrapped `spin_phase`
  accumulator (not a fixed-duration `Animation::Rotation`), so the rate is
  lerpable. `BASE` is Flyby's `make_rotation(Vector(0,0,-1),
  Vector(0,-1,0))`; for wander-dominated presets it is an invisible fixed
  pre-rotation of a random walk.
- **`wander`** — random-walk mix in [0, 1]. Both `RandomWalk`s (inner +
  outer) run permanently on shadow orientations; the effective contribution
  is `slerp(identity, walk_q, wander)`, computed **once per frame**. At 0
  the walks contribute nothing (Flyby); at 1 they contribute fully
  (Liquid2D). Walks keep running while mixed out so a transition into
  wander picks up a live, well-mixed trajectory rather than a cold start.
  Two subtleties are load-bearing:
  - **Hemisphere continuity.** The walk quaternion accumulates unboundedly;
    when its total rotation angle crosses π, `dot(identity, walk_q)` flips
    sign and a naive `slerp(I, walk_q, wander)` jumps between the `w·α` and
    `w·(2π − α)` partial rotations — a visible camera pop for any wander ∈
    (0, 1). Keep a sign-adjusted copy of the walk quaternion (flip its
    stored sign whenever `dot(prev, curr) < 0`) and slerp against that.
  - **Inner mix sense.** Liquid2D applies the inner rotation *forward*
    (`orientation.orient`), while the unified pipeline uses `unorient`
    throughout. An inverted random walk is still a random walk, so this is
    accepted as one of the Liquid2D distribution-level deviations (§10)
    rather than special-cased.

  **Noise generators: three, strictly separated.** `RandomWalk`'s ctor
  reseeds (randomly) and sets the frequency of whatever generator it is
  handed, so each walk owns a private generator and the **warp generator is
  never given to a walk** — Liquid2D's inner walk sharing the warp source
  was an accident that made its warp field randomly seeded and
  frequency-coupled to walk construction order. ShaderBall's warp generator
  gets an explicit `SetFrequency(0.01)` (Flyby's effective default) and the
  default seed, making the warp field deterministic; Liquid2D-mapped preset
  `warp_scale` values are doubled to compensate for the 0.02 → 0.01 base
  frequency (`(z·k)·0.02 ≡ (z·2k)·0.01` bit-for-bit — a ×2 is float-exact).
- **`lens_mix`** — glitch-lens blend in [0, 1], with the shortcuts
  `lens_mix == 0 → rv`, `lens_mix == 1 → glitch(rv)`. The direction-space
  blend `sv = normalize(rv + lens_mix · (glitch(rv) − rv))` specified here is
  **superseded**: `glitch(v) = −v` has exactly two solutions, `(0, 0, ±1)`,
  where the blended direction flips by 180° at `lens_mix = 0.5`, and the
  field is badly distorted across roughly `lens_mix ∈ [0.4, 0.6]`. A
  squared-magnitude guard removes the NaN but not the flip, so the guard is
  not a fix. The shipped blend works in **sample space**: sample the pattern
  once on `rv` and once on `glitch(rv)` and lerp the two `PatternSample`s
  (`blend_lens_samples` in `effects/ShaderBall.h`, which is the source of
  truth for the blend). It is linear in `lens_mix` through the former
  midpoint and costs a second pattern evaluation across the open interval,
  which is the cost §9 budgets.

`lenses::glitch_lens` lives in `core/math/lenses.h` beside its peers. It ships as a
closed-form rational map with no hemisphere fold and no pole guard — its
denominator `1 + 3y²` is at least 1, so the map is finite everywhere on the
sphere and needs neither.

## 6. Color

- **Palette bank.** Bake both existing recipes at init:
  slot 0 = Liquid2D's (STRAIGHT/COMPLEMENTARY/CUP/VIBRANT/75), slot 1 =
  Flyby's (STRAIGHT/SPLIT_COMPLEMENTARY/FLAT/MID/42). A float
  **`palette_pos`** ∈ [0, N−1] selects: integral → single LUT lookup;
  fractional → two lookups + `Color4` lerp. Lerping `palette_pos` between
  presets *is* the palette crossfade — no extra machinery, and the double
  lookup only costs during transitions. The bank is trivially extensible for
  future presets.
- **`breathe_depth`** ∈ [0, 0.3] — `BreatheModifier` amplitude becomes
  preset-driven (0.15 for Liquid2D looks, 0 for Flyby looks). At 0 the
  modifier's coordinate shift is exactly 0, so Flyby lookups are unchanged.
  Driven by `cycle_phase` ← `Driver(cycle_speed)` as in Liquid2D. The
  consuming palette keeps `Wrap = true` (BreatheModifier's `requires_wrap`).
  **Wrap-fold hazard:** `pole_normalize_pattern` can produce exactly 1.0,
  and a Wrap palette folds 1.0 → 0.0 — the palette's *other end* — where
  Flyby's raw `BakedPalette::get` clamps. Clamp `value` to `1.0f − ULP`
  before every lookup (§4) so pattern peaks can't flip color.
- **`hue_shift`** ∈ [0, 1] — Flyby's displacement-driven
  `hue_rotate(−disp · hue_shift)`; skipped when 0 (frame-constant test).
  Note `hue_rotate(c, 0)` is **not** an identity — it round-trips pixel →
  linear RGB → OKLab → gamut map → pixel — and Flyby applies it
  unconditionally, so the skip is a small deliberate visual delta on
  zero-hue-shift presets (F3), accepted as part of the §10 parity softening.
  The gamut LUT is always bought (§8) because any transition may pass
  through nonzero hue_shift.
- **`value_fade`** ∈ [0, 1] — generalizes Flyby's `alpha *= (1 − value)` to
  `alpha *= (1 − value · value_fade)`; 0 reproduces Liquid2D's untouched
  alpha.

## 7. Params, presets, choreography

### Params (all float, all preset-lerped, all registered animated)

| # | Param | Range | Liquid2D value | Flyby value |
|---|---|---|---|---|
| 1 | `warp_scale` | 0.1 – 100 | per preset | per preset |
| 2 | `warp_strength` | 0 – 30 | per preset | per preset |
| 3 | `warp_time_scale` | 0.05 – 1 | 0.5 | 0.3 |
| 4 | `pattern_freq` | 1 – 20 | per preset | per preset |
| 5 | `speed` | 0 – 5 | time_speed | speed |
| 6 | `complexity` | 0 – 3 | 0.5 / 3.0 | 0 |
| 7 | `phase_direct` | 0 – 1 | 0 | 1 |
| 8 | `phase2_rate` | 0 – 2 | 0.8 | 0.7 (was Drift) |
| 9 | `pole_fade` | 1 – 20 | per preset | per preset |
| 10 | `spin_rate` | 0 – 0.05 | 0 | 2π/300 |
| 11 | `wander` | 0 – 1 | 1 | 0 |
| 12 | `lens_mix` | 0 – 1 | 1 | 0 |
| 13 | `palette_pos` | 0 – 1 (N−1) | 0 | 1 |
| 14 | `breathe_depth` | 0 – 0.3 | 0.15 | 0 |
| 15 | `cycle_speed` | 0 – 1 | 0.05 | 0 |
| 16 | `hue_shift` | 0 – 1 | 0 | per preset |
| 17 | `value_fade` | 0 – 1 | 0 | 1 |

Conventions carried over: `params = presets.get()` before any
`register_animated_param` (register captures the default);
`preset_in_ranges` + `all_presets_in_ranges` static_assert; `sizeof(Params)
== 17 * sizeof(float)` static_assert tied to the lerp implementation;
public-first member layout.

Flyby's "Drift" stops being a live-only slider and becomes the
preset-driven `phase2_rate`. That trades a live control for lerpability;
"Pause Animation" + slider takeover still gives live control, matching how
every other preset param already works.

### Preset bank

Seven mapped presets: L0–L1 and F0–F4 carried into the superset via the table
above (unlisted fields take the other family's neutral value), with one mapping
adjustment: Liquid2D `warp_scale` values are **doubled** (1.5 → 3.0) to
compensate for the unified warp generator's 0.01 base frequency (§5).

**Gate-endpoint rule:** the skip branches in §4 stay latched only because
`Animation::Lerp`'s final step lands bit-exactly on the target, which holds
when the gate params' endpoints are exactly 0 or 1 (where `a + (b−a)·t` is
exact) and the easing satisfies `e(1) == 1.0f` (true for `ease_in_out_sin`).
Presets must therefore keep `complexity`, `phase_direct`, `lens_mix`,
`hue_shift`, `value_fade`, and `wander` at exactly 0 or 1 whenever they
intend the corresponding fast path — a near-0 value silently un-latches the
skip forever after the first transition. Pinned by test (§10).

Follow-up presets
should deliberately mix families — e.g. a lensed grid (`lens_mix 0.6,
complexity 0, spin on`) or a hue-warped liquid (`complexity 2, hue_shift
0.4, wander 1`) — that is the payoff of the merge, curated visually after
landing.

### Choreography

A `Choreo` array parallel to `PRESETS` (static_assert equal sizes), one
entry per preset, consumed when *entering* that preset:

```
struct Choreo { uint16_t dwell_min, dwell_max, blend_frames; bool staggered; };
```

Dwell is **sequential**: hold, then blend. The state machine on entering a
preset is:

- `dwell_max == 0`: chain the next `Animation::Lerp` **directly** from the
  previous blend's `.then()` — no timer, zero idle frames (Flyby's exact
  structure). A repeating `RandomTimer(0, 0, …)` is **forbidden** here: it
  fires every step and would spawn a fresh Lerp per frame.
- `dwell_max > 0`: arm a **one-shot** `RandomTimer(dwell_min, dwell_max)`
  whose callback schedules the Lerp; the Lerp's `.then()` re-enters the
  state machine for the next preset.

Entries:

- Flyby-style: `dwell 0/0`, `blend 480`, parallel → perpetual motion.
- Liquid2D-style: `dwell 30/90`, `blend 60`, staggered. (Liquid2D's
  repeating 90–150-frame timer ran the 60-frame blend *inside* its period,
  so the true hold between blends is 30–90 frames; sequential
  dwell-then-blend with 30/90 reproduces the original 90–150-frame cadence.)

Liquid2D's staggered `Params::lerp` (each changed field gets its own time
slice) generalizes unchanged to 17 fields; the `staggered` flag selects
between it and the parallel lerp. Cross-family transitions default to
parallel — a staggered walk through 15 changing fields would spend most of
the blend in unbalanced intermediate states. Scheduling stays on the
existing `timeline.add_pausable` + `Animation::Lerp(...).then(...)` /
`RandomTimer` primitives; both patterns already exist in the two effects.

## 8. State, wrap invariants, memory

Accumulators (all wrapped every frame, exactly the union of today's):
`noise_time` → `STEREO_NOISE_TIME_PERIOD`; `sin_phase`, `phase2`,
`spin_phase`, `cycle_phase` → 2π. The two walk orientations and the
timeline keep Liquid2D's declaration-order rule (orientations before
`timeline`, which clears RandomWalks on teardown).

Persistent-arena footprint:

```
gamut LUT: gamut_lut_bytes(GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS)  // 131,072 B
palettes:  2 × BakedPalette::required_arena_bytes()                   // 2 × 2,052 B
```

`static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET)` as in Flyby
today; the delta vs Flyby is one extra palette bake (+2,052 B), well inside
the partition Flyby already fits.

No heap, arena only; no hardware asserts on the per-pixel path (the nlerp
guard is a branch, not a check).

## 9. Performance budget

Shipping device profiles of the two predecessors (Teensy, 96×20 segmented,
62.5 ms window, 10,368 px quadrant). Both are the last captures taken before
the merge; neither effect is in the roster any more, so they cannot be
regenerated:

- Liquid2D: `lq_shader_draw` 28.3 ms (1,639 cyc/px), ~34 ms headroom
  ([report](../profiles/shipping/profile_liquid2d_teensy_2026-07-25.md)).
- Flyby: `fly_shader_draw` 37.3–44.0 ms across presets, worst preset 44.0
  ms (2,546 cyc/px). The report's headline says 17 ms of margin, but its
  own ISR accounting puts the render budget at ≈59.3 ms/window with the
  peak frame leaving **14.2 ms unspent — that is the binding margin**
  ([report](../profiles/shipping/profile_flyby_teensy_2026-07-27.md)).

Cost deltas for a Flyby-look preset running in ShaderBall (all
frame-constant-branch skips active): one extra orientation transform
(`outer.unorient`: quaternion rotate is 15 mul + 15 add, no division), the
predicted-not-taken lens/complexity/hue branches, and the
`value_fade`/`palette_pos` arithmetic. Estimated well under 100 cyc/px
(< 1.7 ms) against the 14.2 ms margin.

Transient worst case is mid-transition between families: lens nlerp +
cross-coupled pattern (2 extra fast trig) + hue rotate + dual palette
lookup simultaneously. Rough estimate +300–500 cyc/px ≈ +5.2–8.6 ms over
the Flyby worst preset — an order-of-magnitude bound, not a budget, and it
consumes up to ~60% of the honest margin; note the worst Flyby preset is
itself a mid-lerp state, so the stack-up is real. The **first device
profile must capture a cross-family transition window**, not just preset
holds (use the teensy-profile flow; sweep, don't sample).

Mitigations if the transition spills (in order): shorten cross-family
`blend_frames`; schedule cross-family transitions so lens_mix and
complexity ramps don't peak together (stagger just those two fields); only
then consider a second shader variant (ITCM-costly — the ShapeShifter
campaign priced a duplicated inner path at ~7 KB, and the budget has no
room for that by default).

Keep the shader inside an `HS_O3` region exactly as both effects do today
(`stereo_noise_warp` is already `HS_O3_FN`). Re-measure the phantasm size
gate after landing — deleting two effect instantiations and adding one
superset should net out negative, but measure, don't assume.

## 10. Testing and validation

- **Port the white-box wrap test**: seed every accumulator past its ceiling,
  assert `noise_time ∈ [0, STEREO_NOISE_TIME_PERIOD)` and all four trig
  phases ∈ [0, 2π) from frame 0 (union of `test_liquid2d_phase_wrapped` and
  `test_flyby_phase_wrapped`).
- **Port `test_liquid2d_glitch_lens_unit_norm`** unchanged (lens moves
  verbatim), plus a case pinning the nlerp guard fallback.
- **New: formula-reduction test.** Pin `sample(w, φ₁, φ₂)` at
  `(c, s) = (complexity, 0)` against Liquid2D's closed form and `(0, 1)`
  against Flyby's, over a grid of arguments — this is the load-bearing
  equivalence claim of the whole merge. Extend it with the gate-latch
  assertion (§7): simulate a full `Animation::Lerp` between presets and
  assert every gate param lands **bit-exactly** on its 0/1 endpoint, so the
  skip branches provably re-latch after transitions.
- **Visual parity via the framebuffer dump harness** (render → PNG, don't
  theorize). Neither family is bit-exact; do not gate the migration on a
  bit-exact framebuffer diff — it will fail.
  - *Flyby-mapped presets*: near-identical, ULP-level drift. Known exact
    mechanisms: the spin trajectory (Flyby's `Animation::Rotation`
    accumulates incremental quaternion products; ShaderBall rebuilds
    `R_y(spin_phase) ∘ BASE` from a phase accumulator — equal in math,
    different in floats from frame ~2 on); the `hue_rotate` skip at
    `hue_shift == 0` (§6, affects F3); and preset-cycle timing. Compare
    visually / statistically with tight tolerances.
  - *Liquid2D-mapped presets*: distribution-identical only. Deviations:
    walk composition order, the BASE pre-rotation, the inner mix applying
    the walk inverted (`unorient` vs Liquid2D's `orient`, §5), and the warp
    field being deterministic instead of walk-ctor-reseeded (§5).
- Roster oracles update mechanically: smoke suite iterates HS_EFFECT_LIST;
  `GOLDEN_ROSTER` in `tests/test_param_marshal.h` is a hand-maintained
  parallel copy — remove Flyby and Liquid2D rows, insert ShaderBall between
  RingSpin and ShapeShifter.
- Wire new tests into the CI module list and update the roster row in
  `run_tests.cpp` for the deleted Flyby/Liquid2D cases and the added
  ShaderBall ones — the row pins **both** an exact case count (enforced by
  `tests/check_case_calls.cmake`) and the assertion floor, so both columns
  change.

## 11. Roster removal/addition checklist (both repos)

Holosphere:

- `effects/ShaderBall.h` new (`git mv effects/Flyby.h` to keep history on
  the closer parent); Liquid2D's header deleted.
- `core/engine/effects.h`: both `#include`s → one; both `X()` rows → one, at
  the case-insensitive-alphabetical slot (after RingSpin, before
  ShapeShifter) in **both** HS_EFFECT_LIST and HS_PHANTASM_EFFECT_LIST;
  backslashes at column 80. `HS_EFFECT_COUNT` 24→23; the `COUNT − 2`
  phantasm static_assert holds unchanged (both lists shrink by one).
- `tests/test_effects.h`: white-box blocks + runner calls (§10);
  `run_tests.cpp` roster row — case count and assertion floor (§10).
- `tests/test_param_marshal.h`: GOLDEN_ROSTER (§10).
- `README.md`: two gallery cards → one (`?effect=ShaderBall`,
  `docs/screenshots/ShaderBall.png`, alt, heading) **plus** the prose
  references: the pipeline walkthrough naming both effects, the
  `stereo_noise_warp` section, the `hue_rotate` doc, and the gamut-LUT
  passage ("`Flyby` and `MeshFeedback` arm a copy").
- `Holosphere.ino` commented show list.
- `tools/profile_one.sh`: `CYCLERS` list names both effects.
- `tools/profile_sweep.sh`: per-effect sweep rows for both; Flyby's carries
  `-D HS_PROFILE_EPOCH_REVS=2560` — ShaderBall's mixed choreography needs
  its own epoch decision (the preset cycle is longer and partly random).
- `docs/profiles/` READMEs: retire the Liquid2D/Flyby rows and reports; new
  ShaderBall report after the landing profile. Also
  `docs/profiles/memory/arena_high_water.md` (rows for both).
- wasm bootstrap needs no edit (it takes `WASM_EFFECT_NAMES[0]`, which stays
  BZReactionDiffusion).

daydream (separate repo, land via worktree + FF):

- `daydream.js`: Flyby appears only in `HiResFavorites`, Liquid2D only in
  `LoResFavorites` — remove each from its list and add ShaderBall to
  **both** (it spans both looks; the arrays are curated, not alphabetical,
  so the position is a curatorial choice). These hardcoded arrays feed both
  the sim's effect list and the `?effect=` validator; missing this breaks
  `just screenshots`.
- daydream's `README.md` is the Holosphere README mirror — the same gallery
  cards and prose references above must land there too (mirror-accumulator
  flow).
- `daydream/tests/engine_contract_wasm.test.js` needs no edit (verified: it
  names no effect).
- Screenshot order: land rename → install wasm → fix favorites →
  `just screenshots` from the main tree (Playwright lives only there).

## 12. Migration plan

1. **Land ShaderBall alongside both old effects** (roster temporarily 25;
   phantasm assert edited to match for the interim, or do it in one atomic
   land — one atomic land is simpler given the assert). Preferred: single
   commit series in one worktree that adds ShaderBall and removes both, so
   every roster oracle transitions once.
2. Parity validation (framebuffer harness) before the removal commit.
3. Device profile (holds **and** a cross-family transition window) +
   phantasm size gate re-measure.
4. daydream favorites + screenshots.
5. Curate 2–3 new mixed-family presets; retire any that read as redundant.

## 13. Curation decisions

- Bank composition: five wandering-liquid, five spinning-grid, and two
  mixed (grid look on the liquid palette).
- `phase2_rate` ships preset-driven as the animated **Drift** param over
  [0, 2]. No extra live-only slider was added; pause plus slider takeover is
  the live control.
- `spin_rate`'s ceiling ships at the proposed 0.05 rad/frame (≈ 3× Flyby).
