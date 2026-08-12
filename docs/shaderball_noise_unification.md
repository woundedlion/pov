# ShaderBall noise unification

## Decision

Replace ShaderBall's three projected-plane noise warps and the Tangent Noise
lens with one sphere-space `Surface Noise` stage. The stage has `Direct` and
`Curl` motion modes, works with every projection and surface lens, and never
uses a chart coordinate to construct its displacement vector.

Remove Stereo Noise, Vector Noise, and Curl Flow from the Planar Warp 1/2
menus. Keep their numeric enum values as import tombstones, migrate the preset
bank, and delete the legacy stereo kernel after its last compatibility reader
no longer calls it. Exact reproduction of the Stereo Noise look is not a goal.

Move Noise Contour onto the same sphere-field sampler in a later step. It
remains a scalar source, not a displacement mode. Sharing a noise field contract
does not make scalar sampling and vector displacement the same operation.

This design preserves the fixed planar warp program for operations that are
intrinsically planar: Affine Frame, Wave Shear, Vortex, Mirror Tile, and Polar
Chart.

## Goals

- Make noise displacement valid for Folded Sinusoidal, Stereographic,
  Gnomonic, Bonne, Peirce Quincuncial, Airocean, and Equirectangular.
- Remove chart seams, stereographic pole attenuation, and tangent-frame cuts
  from noise displacement.
- Keep noise composable with Glitch, Twist, Mobius, and every kaleidoscope
  lens.
- Give Direct and Curl controls comparable meanings and useful ranges.
- Share basis, seed, animation, channel, resource, and validation contracts
  among displacement noise and Noise Contour.
- Preserve the requested/accepted GUI contract. A rejected edit must not alter
  another requested value or the rendered configuration.
- Fit the Teensy 4.0 device image and frame budget and the WASM segmented-worker
  runtime.

## Non-goals

- Pixel-identical Stereo Noise migration.
- A general node graph or an arbitrary number of noise stages.
- Automatic parameter correction, integrator substitution, or projection
  substitution.
- Making deterministic planar warps projection-independent.
- Replacing FastNoiseLite in the first implementation.

## Current implementations

| Consumer | Coordinate domain | Field construction | Animation | Main limitation |
|---|---|---|---|---|
| Stereo Noise warp | Stereographic plane | Two OpenSimplex2 samples form `dx,dy`; the second channel is offset by 100 on every sampled axis. Generator frequency is 0.01. Strength is multiplied by `1 / (1 + r^2 / pole_fade^2)`. | Dedicated time in `[0,65536]` with a 1024-unit wrap crossfade. UI rate is `[0.05,1]`. | Meaningful only in the stereographic chart. It fades toward the chart's point at infinity, cannot appear twice, and requires Stereographic projection. |
| Vector Noise warp | Selected projection's plane | Two scalar basis samples form a planar vector. The second channel uses fixed offsets. Ridged noise subtracts extra samples to remove its positive DC term. `Vector Angle` rotates the pair. | Stage phase in turns, mapped to a 256-unit noise path with a final 1/64-turn crossfade. | Discontinuous across interrupted or folded chart seams and sensitive to unbounded chart coordinates. |
| Curl Flow warp | Selected projection's plane | A central-difference gradient of one scalar field is rotated 90 degrees. Simplex has a transformed-coordinate fast path. Euler or midpoint integration follows the planar field. | Same stage phase as Vector Noise. | Has the same chart seams. The conservative gradient bound forces `scale * abs(strength) * 64 / intervals <= 0.5`, producing very small useful strengths at ordinary scales. |
| Tangent Noise lens | Unit sphere | Three decorrelated scalar samples form an ambient 3D vector. Its radial component is removed and the result is added to the unit direction and normalized. | Lens phase in turns with the common 256-unit wrap crossfade. | Sphere-safe after removal of the old polar frame switch, but mutually exclusive with every other lens and parameterized differently from the warp noise modes. |
| Noise Contour source | Selected projection's plane | One scalar basis sample, clamped to `[-1,1]`, followed by a contrast remap. | Independent source-noise clock with the common wrap crossfade. | Discontinuous across strict projection seams. It is a value field, not a displacement field. |
| MeshFeedback displacement | Unit sphere | `noise_transform` samples three OpenSimplex2 channels at offsets 0, 100, and 200, multiplies by `amplitude * 0.05`, removes the radial component, caps the tangent chord at 0.5, and renormalizes. | `NoiseParams::time += speed`; time is added to the third sampled coordinate. | The geometric construction is useful, but its scale/frequency split, amplitude calibration, time path, fixed basis, and feedback-specific coarse raster do not match ShaderBall. |

ShaderBall currently stores up to eight prepared `FastNoiseLite` objects. A
`NoiseResourceKey` contains basis, seed, generator-frequency convention,
resource ID, channel version, octave version, and stencil version. Two consumers
with the same resource ID but different keys make the entire requested
configuration invalid.

The existing fixed two-stage planar program evaluates Planar Warp 1 and then
Planar Warp 2 in pullback order. Moving only the noise math while leaving it in those
slots would give the slot names false composition semantics. A sphere
displacement cannot be interleaved with arbitrary plane warps without a valid,
single-valued inverse for every projection.

## Target pipeline

Add one optional sphere-space stage without changing the meaning of the planar
warp positions:

```text
view
  -> inverse outer camera
  -> optional Surface Noise before lens
  -> surface lens
  -> inverse projection frame
  -> optional Surface Noise after lens
  -> sphere-to-plane projection
  -> Planar Warp 1
  -> Planar Warp 2
  -> planar source
```

`Surface Noise Placement` selects which of the two marked positions is active.
There is one field and one displacement application, not two simultaneously
active noise stages.

- `Before Lens` distorts the directions entering the lens. With a
  kaleidoscope, the noise can bend the reflection boundaries and break exact
  repeated-cell symmetry.
- `After Lens` displaces the folded/lensed direction. With a kaleidoscope,
  symmetry-equivalent inputs first collapse to the same chamber direction, so
  the repeated cells retain the same noise detail. This is the recommended
  default.

The placement is a requested setting, never an inferred correction. Selecting
a lens must not change it.

Tangent Noise migrates to `Surface Noise = Direct`. Vector Noise and Curl Flow
migrate to `Direct` and `Curl` respectively. Stereo Noise migrates to `Direct`
and receives preset-specific retuning.

### Lens responsibility

A Lens is a deterministic sphere-domain remap that changes structure,
symmetry, repetition, or topology before projection. Kaleidoscopes fold the
sphere into reflection chambers; Glitch performs latitude doubling and
azimuth tripling; Twist remaps azimuth by latitude; Mobius applies a conformal
sphere transform. Those operations remain lenses because their purpose is the
domain mapping itself.

Stochastic displacement is not a lens. Tangent Noise moves into Surface Noise
even though both stages accept and return a unit direction. The distinction is
semantic and compositional: Lens selects the structural domain, while Surface
Noise perturbs that domain continuously and can be placed explicitly on either
side of the lens.

## Common sphere field

### Field specification

All new consumers use an immutable specification:

```text
NoiseFieldSpec {
  basis
  seed
  scale
  rate
  phase
  channel_layout
  octave_layout
  loop_layout
}
```

`basis` is Simplex, FBM 3, or Ridged 3. Generator frequency is always 1; scale
is a sample input and never mutable generator state. `phase` is runtime state,
not preset state. It advances by the live rate and is forked across discrete
transitions under ShaderBall's existing clock rules.

Use an exactly periodic 3D sampling path:

```text
q(v, phase) = scale * v
            + loop_radius * (cos(2*pi*phase) * A
                           + sin(2*pi*phase) * B)
```

`A` and `B` are fixed orthonormal vectors selected by `loop_layout`.
`loop_radius` and the layout are kernel-version constants, not initial UI
controls. This path is continuous at phase wrap without blending two unrelated
time slices. It also avoids treating one component of the sphere direction as
both spatial Z and time.

FBM 3 uses lacunarity 2, gain 1/2, and normalization by 1.75. Ridged 3 applies
`1 - abs(n)` per octave and removes or recenters its DC term before it is used
as a vector component. These definitions are shared by source and displacement
sampling.

### Direct vector field

Sample three decorrelated channels at `q + C0`, `q + C1`, and `q + C2`:

```text
r = (N(q + C0), N(q + C1), N(q + C2))
u = r - dot(r, v) * v
u = u / max(1, length(u))
u = cos(direction) * u + sin(direction) * cross(v, u)
```

The fixed channel offsets and normalization are part of `channel_layout`.
Projection with `I - vv^T` constructs a tangent vector without choosing a
local east/north frame. It is continuous at both poles and the antimeridian.
The final rotation around `v` gives the current Vector Angle control a
coordinate-free meaning.

### Curl vector field

Curl uses one scalar potential on the sphere. Estimate its ambient gradient
with a symmetric four-sample tetrahedral stencil, project the gradient into the
tangent plane, and rotate it in that plane:

```text
g  = tetrahedral_gradient(N, q, stencil)
gs = g - dot(g, v) * v
u  = cross(v, gs)
u  = u / max(1, length(u))
```

The stencil offsets are applied in noise coordinates. Omitting an additional
factor of `scale` makes Strength an angular travel control instead of allowing
Scale to multiply velocity. The bounded `u` removes the current
scale-strength inequality and its unusably small high-scale range.

The tetrahedral stencil has no pole-dependent basis switch. Its version and
radius are part of the resource key. A central six-sample reference
implementation remains test-only for error comparisons.

### Sphere integration

Move along the sphere with the exponential map rather than chord addition:

```text
exp_v(d) = cos(length(d)) * v + sin(length(d)) * d / length(d)
```

Direct mode applies `exp_v(strength * u(v))` once. Curl mode offers:

| UI value | Field evaluations | Operation |
|---|---:|---|
| Euler | 1 | Evaluate at the input and take one exponential-map step. |
| Midpoint | 2 | Evaluate at the input, evaluate again at the half-step, then take one full step with the midpoint vector. |
| Midpoint 2x | 4 | Perform two midpoint half-steps. |

These names describe field evaluations and substeps directly. Do not retain the
current misleading `Midpoint 2` and `Midpoint 4` labels.

`Strength` is signed for Curl so its sign reverses circulation. Direct uses a
nonnegative Strength and the Direction control for orientation.

## Noise Contour

Noise Contour should ultimately become `Noise Contour (Sphere)`. It samples the
common scalar field at the projection-local sphere direction after the selected
surface lens and Surface Noise stage. The existing contrast transfer remains.

This source intentionally bypasses projection coordinates and the planar warp
program. A requested configuration combining Sphere Noise Contour with a
non-None planar warp is invalid and remains pending with a warning; no warp is
silently disabled. Projection can still supply value weighting and coverage,
but it no longer creates the noise features.

This behavior makes the source continuous for every projection and matches the
fact that one low Scale value covers the whole sphere. The UI label must include
`(Sphere)` so changing Projection and seeing the same underlying noise field is
not surprising.

If chart-space noise remains desirable, add it later as the explicitly named
`Noise Contour (Projected)` source with seam warnings. Do not overload one
option with a hidden domain switch.

## GUI model and ranges

The new controls are:

| Control | Values or range | Notes |
|---|---|---|
| Surface Noise | None, Direct, Curl | Separate from Lens and Planar Warp 1/2. |
| Surface Noise Placement | Before Lens, After Lens | Shown when Surface Noise is active. Default After Lens. |
| Noise Basis | Simplex, FBM 3, Ridged 3 | Shared names and ordering for source and displacement. |
| Noise Scale | `1/64 .. 8` | Logarithmic slider. Scale 1 is approximately sphere-sized structure; the low end is necessary for whole-sphere kaleidoscope looks. |
| Noise Strength | Direct `0 .. 0.5` rad; Curl `-0.5 .. 0.5` rad | At most about 29 degrees of surface travel per application. |
| Noise Rate | `-1/1024 .. 1/1024` turns/frame | Fine primary slider. Exact entry may accept up to `+/-1/64` when performance and alias tests pass. |
| Noise Direction | `0 .. 1` turn | Direct only; rotates the tangent vector around the sample direction. |
| Curl Integrator | Euler, Midpoint, Midpoint 2x | Curl only. No automatic quality change. |
| Noise Seed | Signed 32-bit integer | Exported and preset-authored; expose as exact entry only if the GUI supports integer editing. |

There is no Pole Fade, projection envelope, or Edge Fade in the sphere-noise
stage. Those are chart policies. Projection coverage and value weighting remain
independent downstream controls.

Keep source labels prefixed with `Source` and surface-displacement labels
prefixed with `Surface`. Dynamic schemas must not reuse one visible name for two
different targets.

The GUI presents one collapsible bank per pipeline responsibility. In display
order: Camera, Surface Noise, Lens, Projection Frame, Projection, Planar Warp
1, Planar Warp 2, Function, Signal Weight, Value Transfer, Coverage, and Color.
Surface Noise Placement states whether that bank executes immediately before
or after Lens; the UI never moves another bank or changes placement as a side
effect. Bank organization is URL-transparent, and renamed legacy keys migrate
to their canonical names without changing values.

## Requested and accepted state

Noise unification retains ShaderBall's two-state admission model:

- `requested_config` is exactly what the user selected, including invalid
  combinations.
- `accepted_config` is the last complete valid configuration and is the only
  configuration that may reach frame preparation and rendering.
- A single edit writes only its requested target. It must not alter Projection,
  Lens, Surface Noise Placement, integrator, scale, strength, either planar
  warp, or any hidden accepted value.
- If several pending edits make the complete requested configuration valid,
  admit that configuration atomically.
- Rejection leaves rendering byte-for-byte on the prior accepted
  configuration. There is no fallback noise kind.
- The GUI displays requested values and warning markers. Renderer snapshots,
  reload seeds, and worker rebuilds use accepted values first, then replay
  requested values.

Every warning names the failing controls, prints their requested values, states
the exact range or relationship, and gives actions that can make the whole
request valid. Examples:

```text
Surface Noise Curl rejected: Noise Strength 0.7 is outside
[-0.5, 0.5]. Set Noise Strength within that range.
```

```text
Noise Contour (Sphere) cannot be combined with Planar Warp 1 Mirror Tile.
Set Planar Warp 1 to None or choose a planar Function.
```

Warnings must not recommend an automatic integrator upgrade. Selecting a
higher-quality integrator is an action the user may take, not a mutation the
engine may perform.

Browser tests must exercise the real control event path. Native admission tests
alone cannot detect a stale HTML select after the engine has accepted a value.
At minimum, the browser suite must select Lens Glitch and then Lens None and
assert that the rendered engine value, selected option, URL value, and worker
snapshot all say None.

## Resource ownership

Replace consumer-specific resource keys with one `NoiseFieldKey` containing:

- basis and seed;
- scalar, direct-channel, or curl-stencil layout version;
- octave normalization version;
- periodic loop layout and radius version; and
- any FastNoiseLite type or generator-frequency convention.

Scale, rate, phase, strength, direction, and integrator are sample/runtime
inputs and do not belong in the generator key. No pixel path may call a
FastNoiseLite setter.

Remove author-selected numeric resource IDs from new presets. Deduplicate
prepared owners by the full key. Resource identity is an implementation detail,
not a visual control. Preserve old resource IDs only while reading legacy
configs.

With one Surface Noise stage and one sphere Noise Contour source, one endpoint
needs at most two field owners. A through-clear transition needs at most four,
well below the existing fixed capacity of eight. Keep the capacity until an ELF
and RAM measurement proves reducing it is useful.

Mutable phase belongs to each consumer clock, not to the shared generator. Two
consumers may share a prepared generator while advancing different phases.

## Presets, enums, exports, and URLs

Numeric compatibility and visual compatibility are separate requirements.

1. Add a serialized ShaderBall schema version before changing option meaning.
2. Add a new Surface Noise enum and parameter block. Do not reuse the numeric
   value of Tangent Noise or a planar warp kind.
3. Keep `LEGACY_STEREO_NOISE`, `VECTOR_NOISE`, `CURL_FLOW`, and
   `TANGENT_NOISE` numeric values reserved as import-only tombstones. They do
   not appear in GUI option arrays or new exports.
4. Migrate every built-in preset by hand and review its rendered result. Stereo
   Noise presets may become Direct, Curl, or no noise; preserving variety is
   more important than mechanical one-to-one mapping.
5. Import old configs into the new model before admission. Migrate accepted and
   requested snapshots independently. An invalid old requested value must not
   overwrite its valid accepted snapshot.
6. Show one non-blocking import notice that names each transformed legacy
   option. This is serialization migration, not a live control side effect.
7. Emit only the new schema from Export and URL writers. Preserve the old
   reader and enum tombstones as long as stored links are supported.

Suggested initial migration heuristics are only starting points for preset
retuning:

- Tangent Noise -> Surface Noise Direct, After Lens; keep basis, seed, scale,
  rate, direction 0, and clamp Amount into the new Strength range.
- Vector Noise -> Surface Noise Direct, After Lens; keep basis, seed, scale,
  rate, and Vector Angle; calibrate Strength in radians from reference renders.
- Curl Flow -> Surface Noise Curl, After Lens; keep basis, seed, scale, signed
  strength, and rate; map Euler 1 to Euler, and both old midpoint settings to
  the nearest measured new quality.
- Stereo Noise -> Surface Noise Direct, After Lens; begin with
  `new_scale = clamp(old_scale * 0.01, 1/64, 8)` and retune Strength per preset.
  Pole Fade is not migrated into the noise stage.

Do not silently map an obsolete live menu selection while the application is
running. Migration occurs only while decoding a versioned stored config.

## Teensy and WASM constraints

ShaderBall runs per visible sphere sample at 288x144 on a 600 MHz Cortex-M7.
The current shipping roster has only hundreds of bytes of ITCM headroom, and
the existing profile already contains red transition windows. New templates or
always-inline specializations cannot be assumed to fit.

Approximate raw noise sample counts per pixel are:

| Mode | Simplex | FBM 3 |
|---|---:|---:|
| Direct | 3 | 9 |
| Curl Euler | 4 | 12 |
| Curl Midpoint | 8 | 24 |
| Curl Midpoint 2x | 16 | 48 |

These counts exclude projection, lens, and source sampling. Ridged 3 has FBM 3
octave count and may require additional recentering work. Midpoint 2x with FBM
3 must not ship on Teensy solely because it is acceptable in WASM.

The first implementation should keep one dispatching kernel in flash and reuse
the existing selectively optimized OpenSimplex2 primitive. Measure before
adding specialization. If direct per-pixel Curl misses the device gate, build a
coarse sphere field once per frame using the existing latitude-ring
`SphericalField` layout:

- store ambient 3D vectors as snorm16 triples;
- interpolate in the ring field with longitude wrap;
- reproject the interpolated vector onto the query point's tangent plane; and
- keep pole samples single-valued.

This follows MeshFeedback's successful coarse-field architecture without
sharing its feedback-specific coordinate projection or parameter semantics.
Allocate the field from a bounded scratch scope. Do not add an unmeasured
persistent array for every possible transition owner.

Device admission may reject a basis/integrator pair whose measured worst-case
cost exceeds the target, but the warning must name that pair and suggest an
explicit cheaper choice. It must not change the pair. Prefer a build-time
capability table over a runtime timing guess. WASM may expose additional
quality only if exports encode it and device imports report the unsupported
choice accurately.

WASM kernels must remain allocation-free per pixel and deterministic with the
native build. Segmented workers receive accepted and requested values
separately under the current protocol. A worker rebuild restores the complete
accepted snapshot before replaying pending requests.

## Implementation stages

### Stage 1: shared field primitives

- Add `NoiseFieldSpec`, `NoiseFieldKey`, periodic 3D sampling, common octave
  normalization, direct tangent construction, and the tetrahedral gradient.
- Add mathematical host tests against scalar reference implementations.
- Do not expose UI options or migrate presets yet.

### Stage 2: Surface Noise

- Add the dedicated stage, placement, Direct mode, Curl mode, exponential-map
  stepping, and explicit integrator options.
- Make it composable with every lens and projection.
- Add native, WASM-contract, worker, and real browser-control tests.
- Profile each basis/integrator/placement family on device before adding it to
  a preset.

### Stage 3: preset and serialization migration

- Add the schema version and legacy readers.
- Hand-retune the preset bank, emphasizing new looks rather than parity.
- Remove obsolete options from GUI arrays and new exports.
- Retain enum tombstones and import tests.

### Stage 4: Noise Contour sphere migration

- Add the projection-local sphere direction to the source sample input.
- Rename the source to `Noise Contour (Sphere)` and use the common scalar
  sampler.
- Reject non-None planar warps without changing them.
- Retune low-frequency kaleidoscope presets.

### Stage 5: legacy deletion and optimization

- Delete `stereo_noise_warp`, the legacy time clock, pole-coupled warp
  validation, and legacy-only prepared-state fields after no runtime path calls
  them.
- Delete projected Vector/Curl kernels after import migration is config-only.
- Remove Tangent Noise lens runtime code after its tombstone is import-only.
- Re-run full roster ELF, arena, native, WASM, and on-device profile gates.

## Tests and acceptance criteria

### Geometry

- Direct and Curl output unit vectors within `1e-5` at both poles, axes, and a
  dense deterministic sphere sample.
- Samples approaching a pole or antimeridian from either side converge within
  a tolerance derived from the noise stencil.
- Angular displacement never exceeds `abs(Strength)` plus float tolerance.
- Period phase 0 equals phase 1 for every basis and channel layout.
- Tetrahedral Curl agrees with the six-sample reference within a measured
  angular-error bound and has no persistent radial component.
- Negative Curl Strength reverses local circulation.

### Composition

- Every Surface Noise mode, basis, placement, and admitted integrator produces
  finite output with all seven projections.
- The same matrix passes with None, Glitch, Twist, Mobius, and all
  kaleidoscope solids represented in the current lens catalog.
- `After Lens` preserves equality for symmetry-equivalent kaleidoscope inputs.
- Interrupted projection layouts do not reject Surface Noise and do not add a
  noise-value discontinuity at glued edges.
- Projection-specific cut coverage may still be discontinuous where the
  projection intentionally cuts the image; that is not a noise seam.

### State and GUI

- Pairwise selector tests cover every Function, Projection, Lens, Surface
  Noise, placement, Planar Warp 1, and Planar Warp 2 choice.
- Changing any selector changes exactly one requested field.
- Rejected requests leave every accepted field and rendered frame state
  unchanged.
- A jointly repaired requested configuration admits atomically.
- Lens Glitch -> None updates the real dropdown, engine requested value,
  accepted value, URL, worker snapshot, and render.
- Reload and worker reconstruction preserve an accepted None plus a rejected
  requested noise value without falling back to a preset.
- Every rejection path has an exact, actionable warning and no generic
  fallback in expected configurations.

### Compatibility

- Every legacy enum number decodes deterministically.
- Accepted and requested legacy snapshots migrate independently.
- New export -> import is lossless for all Surface Noise parameters.
- Old imports produce a visible migration notice and never emit an obsolete
  option on the next export.
- All migrated presets are valid holds and contain no import tombstones.

### Performance and memory

- Native ShaderBall and parameter-marshal suites pass.
- WASM smoke passes at 96x20 and 288x144 for all effects under the CSP and stack
  gates.
- Daydream unit, browser-control, worker, typecheck, and lint suites pass.
- The full Phantasm ELF remains below FLASH, ITCM, and RAM1 region limits.
- Arena high-water remains below the device partition limits.
- Fixed-preset and transition device captures introduce no new red preset
  family. Any existing red family must not regress by more than the profiling
  tolerance agreed for that capture.
- Stage counters report field construction, interpolation, Direct, and each
  Curl integrator separately so a regression can be attributed.

## Risks and open choices

### Curl quality on Teensy

Midpoint 2x with a three-octave basis may be too expensive per pixel. The
recommended order is direct measurement, then a shared coarse sphere field,
then a device capability restriction. Removing Curl or silently downgrading it
is not recommended.

### Sphere Noise placement

One placement selector provides two materially different lens interactions
without the memory and UI cost of two simultaneous sphere stages. If authored
looks later demonstrate a need for two, add a second explicit stage only after
measuring its transition resource union and frame cost.

### Noise Contour and planar warps

A sphere scalar source cannot honestly participate in the current planar warp
stack. The recommendation is to reject that combination and keep a separately
named projected source only if artists still need it. Ignoring the planar warps
would violate the no-side-effect contract at the rendered-behavior level.

### Temporal motion

The circular 3D domain path is exactly periodic and inexpensive, but it is not
4D simplex evolution. Adopt it first. A future 4D kernel requires separate
quality, cost, resource-key, and migration review.

### Migration strength

There is no global conversion from planar stereo displacement to angular
sphere displacement because the old angular effect depends on chart radius and
Pole Fade. Hand-retuned presets are the source of truth. The importer heuristic
only prevents old links from failing to load.

## Recommendation summary

Ship one dedicated sphere-space Surface Noise stage with Direct and Curl modes,
an explicit before/after-lens placement, angular Strength, low whole-sphere
Scale, and bounded coordinate-free vector fields. Remove all noise choices from
the planar warp menus and migrate Tangent Noise out of the Lens menu. Move Noise
Contour to the shared sphere scalar field after displacement migration is
stable. Preserve old numeric values only as versioned import tombstones.

This produces more combinations than Stereo Noise while eliminating its
projection dependency. It also makes invalid state handling simpler: normal
noise choices no longer need projection, pole, seam, coordinate-bound, or
scale-strength rejection rules.
