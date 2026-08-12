# ShaderBall noise unification

## Decision

Replace legacy Stereo Noise and the Tangent Noise lens with one sphere-space
`Surface Noise` stage. The stage has `Direct` and `Curl` motion modes, works
with every projection and surface lens, and never uses a chart coordinate to
construct its displacement vector.

Keep Vector Noise and Curl Flow as explicitly named `Projected Vector Noise`
and `Projected Curl Flow` choices in both planar-warp menus. They remain useful
chart-space effects and may be combined with each other or other planar warps.
Remove only Stereo Noise from those menus, reserve its numeric enum value as an
import tombstone, migrate the preset bank, and delete the legacy stereo kernel
after its last compatibility reader no longer calls it. Exact reproduction of
the Stereo Noise look is not a goal.

Keep the existing source as `Noise Contour (Projected)` with its serialized
value and behavior unchanged. Add `Noise Contour (Sphere)` on the shared
sphere-field sampler as a new source value. Both remain scalar sources, not
displacement modes. Sharing a noise-field contract does not make scalar
sampling and vector displacement the same operation.

This design preserves the fixed planar warp program for operations that are
intrinsically planar: Affine Frame, Wave Shear, Vortex, Projected Vector Noise,
Projected Curl Flow, Mirror Tile, and Polar Chart.

Planar Warp 1 and Planar Warp 2 remain independent authored stages. Every
compatible pair runs sequentially and exposes separate parameters; adding
sphere-space noise does not collapse the planar program. Polar Chart remains
the explicit sole-stage exception because it changes the coordinate
interpretation consumed by the source.

## Goals

- Make sphere-space noise displacement valid for Folded Sinusoidal, Stereographic,
  Gnomonic, Bonne, Peirce Quincuncial, Airocean, and Equirectangular.
- Provide noise displacement without chart seams, stereographic pole
  attenuation, or tangent-frame cuts.
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
- Making projected planar warps projection-independent.
- Replacing FastNoiseLite in the first implementation.

## Current implementations

| Consumer | Coordinate domain | Field construction | Animation | Main limitation |
|---|---|---|---|---|
| Stereo Noise warp | Stereographic plane | Two OpenSimplex2 samples form `dx,dy`; the second channel is offset by 100 on every sampled axis. Generator frequency is 0.01. Strength is multiplied by `1 / (1 + r^2 / pole_fade^2)`. | Dedicated time in `[0,65536]` with a 1024-unit wrap crossfade. UI rate is `[0.05,1]`. | Meaningful only in the stereographic chart. It fades toward the chart's point at infinity, cannot appear twice, and requires Stereographic projection. |
| Vector Noise warp | Selected projection's plane | Two scalar basis samples form a planar vector. The second channel uses fixed offsets. Ridged noise subtracts extra samples to remove its positive DC term. `Vector Angle` rotates the pair. | Stage phase in turns, mapped to a 256-unit noise path with a final 1/64-turn crossfade. | Discontinuous across projection cut boundaries or folded chart seams and sensitive to unbounded chart coordinates. |
| Curl Flow warp | Selected projection's plane | A central-difference gradient of one scalar field is rotated 90 degrees. Simplex has a transformed-coordinate fast path. Euler or midpoint integration follows the planar field. | Same stage phase as Vector Noise. | Has the same chart seams. The conservative gradient bound forces `scale * abs(strength) * 64 / intervals <= 0.5`, producing very small useful strengths at ordinary scales. |
| Tangent Noise lens | Unit sphere | Three decorrelated scalar samples form an ambient 3D vector. Its radial component is removed and the result is added to the unit direction and normalized. | Lens phase in turns with the common 256-unit wrap crossfade. | Sphere-safe after removal of the old polar frame switch, but mutually exclusive with every other lens and parameterized differently from the warp noise modes. |
| Noise Contour (Projected) source | Selected projection's plane | One scalar basis sample, clamped to `[-1,1]`, followed by a contrast remap. | Independent source-noise clock with the common wrap crossfade. | Discontinuous across cut-layout projection seams. It is a value field, not a displacement field. |
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

Tangent Noise migrates to `Surface Noise = Direct`. Stereo Noise migrates to
`Direct` and receives preset-specific retuning. Projected Vector Noise and
Projected Curl Flow remain in their original planar slots and retain their
coordinate-space behavior.

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

Remove the universal Lens Mix control. A structural mapping is either `None`
or fully active; a fractional kaleidoscope or Glitch mapping has no stable
topological meaning. Lens selection and Surface Noise placement changes use
ShaderBall's through-clear transition: render only the old endpoint while it
fades to the clear frame, then render only the new endpoint while it fades in.
Continuously parameterized lenses expose controls with lens-specific meaning.
No serialized or live configuration contains a blended sphere direction
between identity and a lens, and no frame shades both structural endpoints.

## Common sphere field

### Field specification

All noise consumers use an immutable specification:

```text
NoiseFieldSpec {
  domain
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

`domain` is `Sphere3D` or `Projected2D` and is part of the resource key.
Surface Noise and Noise Contour (Sphere) use the sphere construction below.
Projected Vector Noise and Projected Curl Flow retain their two-dimensional
vector and gradient geometry, including chart-boundary warnings, while sharing
the basis, seed, octave, periodic animation, key, and prepared-owner contract.

`basis` is Simplex, FBM 3, or Ridged 3. Generator frequency is always 1; scale
is a sample input and never mutable generator state. `phase` is runtime state,
not preset state. It advances by the live rate and is forked across discrete
transitions under ShaderBall's existing clock rules.

Version 1 uses an exactly periodic three-dimensional sampling path. First set
`t = wrap_t(phase)` and `phi = 2*pi*t`, then construct the noise coordinate for
the consumer's declared domain:

```text
q_sphere(v, phase) = scale * v
                    + 32 * (cos(phi) * (1,0,0)
                          + sin(phi) * (0,1,0))

q_projected(p, phase) = (scale*p.x, scale*p.y, 0)
                       + 32 * (cos(phi) * (0,0,1)
                             + sin(phi) * (1,1,0)/sqrt(2))
```

These planes and radius are the exact `loop_layout = 1` contract for their
domains, not initial UI controls. `wrap_t` is applied before either trigonometric
call, so phase 0 and every integer phase use the identical input coordinate.
This path is continuous at phase wrap without blending two unrelated time
slices. It also avoids treating one component of the sphere direction as both
spatial Z and time.

Version 1 channel offsets are `C0 = (0,0,0)`, `C1 = (17,-29,43)`, and
`C2 = (-47,11,-23)`. Sphere Direct uses all three; Projected Vector Noise uses
`C0` for X and `C1` for Y; scalar and Curl consumers use `C0`. The values,
channel ordering, and Ridged recentering are part of `channel_layout = 1` and
must be golden-tested rather than inferred from consumer names.

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
theta = 2*pi*direction
u = cos(theta) * u + sin(theta) * cross(v, u)
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

The version 1 tetrahedral stencil uses radius `h = 1/64`, directions
`(1,1,1)`, `(1,-1,-1)`, `(-1,1,-1)`, and `(-1,-1,1)`, each divided by
`sqrt(3)`, and computes `g = 3/(4*h) * sum(N(q+h*Di)*Di)`. It has no
pole-dependent basis switch. Its version and radius are part of the resource
key. A central six-sample reference implementation remains test-only for error
comparisons.

Projected Curl Flow differentiates only the two planar coordinates while
holding phase fixed:

```text
nx = (N(q_projected + h*(1,0,0)) - N(q_projected - h*(1,0,0))) / (2*h)
ny = (N(q_projected + h*(0,1,0)) - N(q_projected - h*(0,1,0))) / (2*h)
u  = (ny, -nx)
```

It uses the same version 1 `h = 1/64`. Vector Angle rotates a Projected Vector
Noise pair in the plane after channel sampling. These definitions close the
periodic-field contract for both coordinate domains; no consumer supplies an
implicit time axis or private wrap crossfade.

### Sphere integration

Move along the sphere with the exponential map rather than chord addition:

```text
exp_v(d) = cos(length(d)) * v + sin(length(d)) * d / length(d)
```

Define `exp_v(0) = v`. The implementation branches below a fixed epsilon or
uses a stable `sinc(length(d))` evaluation, so zero Strength and a zero field
are exact identity operations rather than `0/0`.

Direct mode applies `exp_v(strength * u(v))` once. Curl mode offers:

| UI value | Field evaluations | Operation |
|---|---:|---|
| Euler | 1 | Evaluate at the input and take one exponential-map step. |
| Midpoint | 2 | Evaluate at the input, step to the geodesic midpoint, evaluate there, parallel-transport that tangent vector back to the substep start, then take the full exponential-map step. |
| Midpoint 2x | 4 | Perform two midpoint substeps of half the total signed Strength, transporting within each substep. |

For a substep from `v` through midpoint `m`, transport the midpoint tangent
`w` back along the same shortest geodesic before applying it at `v`:

```text
transport_m_to_v(w) = w - dot(w, v) / (1 + dot(m, v)) * (m + v)
```

The maximum substep is below the antipodal singularity. Every vector passed to
`exp_v` belongs to `T_v S^2`; a vector evaluated in `T_m S^2` is never applied
directly at `v`.

These names describe field evaluations and substeps directly. Do not retain the
current misleading `Midpoint 2` and `Midpoint 4` labels.

`Strength` is signed for Curl so its sign reverses circulation. Direct uses a
nonnegative Strength and the Direction control for orientation.

## Noise Contour

Keep the current source, numeric value, and sampling behavior under the explicit
label `Noise Contour (Projected)`. Add `Noise Contour (Sphere)` as a new enum
value. The new source samples the common scalar field at the sphere direction
after the selected surface lens and Surface Noise stage. Both retain the
existing contrast transfer.

This source intentionally bypasses projection coordinates and the planar warp
program. A requested configuration combining Sphere Noise Contour with a
non-None planar warp is invalid and remains pending with a warning; no warp is
silently disabled. Projection can still supply value weighting and coverage,
but it no longer creates the noise features.

Sphere behavior makes the source continuous for every projection and matches the
fact that one low Scale value covers the whole sphere. The UI label must include
`(Sphere)` so changing Projection and seeing the same underlying noise field is
not surprising.

Projected behavior remains useful when planar warps should shape the scalar
field. It follows the same projection compatibility matrix as Projected Vector
Noise and Projected Curl Flow. Do not overload one option with a hidden domain
switch, migrate an existing Projected selection to Sphere, or silently bypass
either planar warp.

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

`getParameterDefinitions()` is a dynamic presentation schema and is not a
persistence format. Add a versioned full-config snapshot API that serializes
all accepted and requested fields, including controls inactive in the current
schema, plus stable field identifiers for pending edits. Worker initialization,
reload, URL persistence, and Export use this snapshot. Restore validates and
installs the complete accepted config atomically, installs the requested config
without correction, reconstructs pending edits by field identifier, and only
then rebuilds the dynamic GUI schema. Per-control replay order is not part of
the restore contract.

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

### Rejection boundary

Surface Noise has no cross-stage Lens or Projection rejection. Direct, Curl,
and both placements compose with every lens and projection. Its admission
rules are limited to its own parameter ranges and explicitly published device
capabilities such as a basis-integrator pair.

Cross-stage rejection remains only where coordinate semantics require it:

- projected-noise consumers follow the exact projection matrix below;
- Polar Chart requires a compatible periodic Function and must be the only
  active planar warp; and
- Noise Contour (Sphere) rejects non-None planar warps because a sphere scalar
  source cannot consume coordinates that exist only in a projection plane.

Signal Weight, Coverage, Color, Lens, and Surface Noise do not acquire hidden
compatibility dependencies. The admitted maximum resource union fits the fixed
capacity, so ordinary stage combinations never surface a resource-capacity
warning.

### Projected-noise projection matrix

The following rule is normative for each Projected Vector Noise or Projected
Curl Flow stage and for the Noise Contour (Projected) source. Layout or
hemisphere suboptions do not alter the result.

| Projection | Projected-noise consumer | Reason |
|---|---|---|
| Folded Sinusoidal | Admit | The fold and its coverage are explicit authored chart behavior. |
| Stereographic | Admit | One unbounded chart with an explicit pole policy. |
| Gnomonic: Folded, Front, or Back | Admit | The selected hemisphere policy explicitly owns the cut, clipping, and singular boundary. |
| Bonne: North or South | Reject | Its cut-layout chart is not a single connected plane for a continuous projected field. |
| Peirce Quincuncial: every layout | Reject | Its square/diamond edge identifications are not represented in the planar sampler. |
| Airocean: Vertical or Horizontal | Reject | Disconnected face layout boundaries are not represented in the planar sampler. |
| Equirectangular | Admit | Its longitude cut is explicit and coordinates remain bounded. |

An invalid warning identifies every offending stage or source. For example:

```text
Planar Warp 2 Projected Curl Flow cannot use Peirce Quincuncial (Diamond):
the projected field does not join across that cut layout. Choose Folded
Sinusoidal, Stereographic, Gnomonic, or Equirectangular; or change Planar Warp 2.
```

The engine never substitutes a projection, a planar warp, or the Sphere source.
Independent per-mode coordinate and parameter bounds may still reject a value
after the matrix admits the projection; that warning must name the bound and
the exact control that can satisfy it.

Browser tests must exercise the real control event path. Native admission tests
alone cannot detect a stale HTML select after the engine has accepted a value.
At minimum, the browser suite must select Lens Glitch and then Lens None and
assert that the rendered engine value, selected option, URL value, and worker
snapshot all say None.

## Resource ownership

Replace consumer-specific resource keys with one `NoiseFieldKey` containing:

- basis and seed;
- sphere-3D or projected-2D coordinate-domain layout;
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

With Surface Noise, both projected planar warps, and Noise Contour (Projected),
one endpoint needs at most four field owners. A through-clear transition may
prepare the union of both endpoints and therefore needs at most eight, exactly
the existing fixed capacity, but the renderer samples at most one four-owner
endpoint in a frame. Admission tests must exercise the eight-distinct-key
preparation maximum and reject only a malformed or future schema that actually
exceeds it. Keep the capacity until an ELF and RAM measurement proves changing
it is useful.

Mutable phase belongs to each consumer clock, not to the shared generator. Two
consumers may share a prepared generator while advancing different phases.

## Presets, enums, exports, and URLs

Numeric compatibility and visual compatibility are separate requirements.

Serialized enum storage IDs are not GUI ordinals. Each visible enum option is
defined as `{storage_id, label, export_name}`. The parameter bridge presents a
dense zero-based list to the dropdown and maps its visible ordinal to the
option's storage ID before writing requested state; the reverse mapping selects
the GUI option from requested state. Snapshots, URLs, presets, and exports store
only storage IDs. No code casts a visible ordinal directly to an enum after an
import tombstone creates a hole.

An import-only tombstone participates only in the versioned decoder. It is
absent from the dense visible list, cannot be selected or emitted, and produces
the documented migration result or an explicit unsupported-value error. Tests
must select every visible option on both sides of each numeric hole and assert
the requested storage ID, displayed ordinal, export name, URL, and restored
selection.

1. Add a serialized ShaderBall schema version before changing option meaning.
2. Add a new Surface Noise enum and parameter block. Do not reuse the numeric
   value of Tangent Noise or a planar warp kind.
3. Keep `LEGACY_STEREO_NOISE` and `TANGENT_NOISE` numeric values reserved as
   import-only tombstones. `VECTOR_NOISE` and `CURL_FLOW` retain their numeric
   planar-warp values and export as Projected Vector Noise and Projected Curl
   Flow.
4. Keep the existing Noise Contour storage ID as Noise Contour (Projected), and
   append Noise Contour (Sphere) with a new ID. Existing snapshots therefore do
   not change noise domain during decoding.
5. Migrate every built-in preset by hand and review its rendered result. Stereo
   Noise presets may become Direct, Curl, or no noise; preserving variety is
   more important than mechanical one-to-one mapping.
6. Import old configs into the new model before admission. Migrate accepted and
   requested snapshots independently. An invalid old requested value must not
   overwrite its valid accepted snapshot.
7. Show one non-blocking import notice that names each transformed legacy
   option. This is serialization migration, not a live control side effect.
8. Emit only the new schema from Export and URL writers. Preserve the old
   reader and enum tombstones as long as stored links are supported.
9. Remove Lens Mix from the new schema. A legacy Lens Mix of zero imports as
   Lens None. A nonzero mix imports a selected structural lens at full effect.
   Tangent Noise with zero mix has no Surface Noise effect; with nonzero mix it
   is a migration candidate whose initial Strength is calibrated from Amount
   and Mix before hand-retuning. Built-in presets are hand-retuned, and external
   imports receive a notice that the fractional blend was removed.

A legacy snapshot may contain Tangent Noise and one Stereo Noise stage
simultaneously, while the new model contains one Surface Noise stage. Migrate
each accepted and requested snapshot independently. An effective Tangent Noise
lens has precedence over Stereo Noise because it already operates on the
sphere. Otherwise Stereo Noise becomes Surface Noise. Clear only the migrated
Stereo Noise slot and the Tangent Noise lens selection; Projected Vector Noise,
Projected Curl Flow, and every non-noise planar stage remain in place. The
import notice names the selected source, a discarded Stereo Noise candidate,
each cleared obsolete slot, and the Lens Mix conversion. This deliberately
lossy collision policy is version-decoding behavior, never a live menu side
effect, and does not promise to reproduce the old look.

Suggested initial migration heuristics are only starting points for preset
retuning:

- Tangent Noise -> Surface Noise Direct, After Lens; keep basis, seed, scale,
  rate, direction 0, and clamp Amount into the new Strength range.
- Vector Noise -> Projected Vector Noise in the same planar slot with unchanged
  parameters.
- Curl Flow -> Projected Curl Flow in the same planar slot with unchanged
  parameters and integrator.
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
native build. Replace the worker protocol's dynamic parameter list with the
versioned full-config snapshot. A worker rebuild restores the complete accepted
config atomically, then restores requested state and pending field identifiers.

## Implementation stages

### Stage 1: shared field primitives

- Add `NoiseFieldSpec`, `NoiseFieldKey`, periodic 3D sampling, common octave
  normalization, direct tangent construction, and the tetrahedral gradient.
- Add mathematical host tests against scalar reference implementations.
- Do not expose UI options or migrate presets yet.

### Stage 2: Surface Noise

- Add the dedicated stage, placement, Direct mode, Curl mode, exponential-map
  stepping, and explicit integrator options.
- Remove Lens Mix from the authored schema and route lens changes through
  through-clear endpoint transitions.
- Make it composable with every lens and projection.
- Add native, WASM-contract, worker, and real browser-control tests.
- Profile each basis/integrator/placement family on device before adding it to
  a preset.

### Stage 3: preset and serialization migration

- Add the schema version, full accepted/requested config snapshot API, stable
  pending-field identifiers, worker protocol, and legacy readers.
- Hand-retune the preset bank, emphasizing new looks rather than parity.
- Remove obsolete options from GUI arrays and new exports.
- Retain enum tombstones and import tests.

### Stage 4: Noise Contour sphere option

- Add the projection-local sphere direction to the source sample input.
- Rename the existing source to `Noise Contour (Projected)` without changing
  its storage ID or behavior.
- Append `Noise Contour (Sphere)` and use the common scalar sampler.
- Reject non-None planar warps for the Sphere option without changing them;
  apply the projected-noise matrix to the Projected option.
- Retune low-frequency kaleidoscope presets.

### Stage 5: legacy deletion and optimization

- Delete `stereo_noise_warp`, the legacy time clock, pole-coupled warp
  validation, and legacy-only prepared-state fields after no runtime path calls
  them.
- Retain the projected Vector/Curl kernels under their explicit projected
  labels and move their basis, key, and animation plumbing onto the shared
  noise infrastructure.
- Remove Tangent Noise lens runtime code after its tombstone is import-only.
- Re-run full roster ELF, arena, native, WASM, and on-device profile gates.

## Tests and acceptance criteria

### Geometry

- Direct and Curl output unit vectors within `1e-5` at both poles, axes, and a
  dense deterministic sphere sample.
- Samples approaching a pole or antimeridian from either side converge within
  a tolerance derived from the noise stencil.
- Angular displacement never exceeds `abs(Strength)` plus float tolerance.
- Direction 0 and Direction 1 produce identical output, and zero Strength is
  an exact finite identity for every basis and integrator.
- Each Euler and midpoint step passes a vector tangent to the step's own base
  point into the exponential map; all integrators preserve unit length.
- Period phase 0 equals phase 1 bit-for-bit after coordinate construction for
  Sphere3D and Projected2D, every basis, scalar/direct/vector/curl channel
  layout, and both projected Curl derivative axes.
- Tetrahedral Curl agrees with the six-sample reference within a measured
  angular-error bound and has no persistent radial component.
- Negative Curl Strength reverses local circulation.

### Composition

- Every Surface Noise mode, basis, placement, and admitted integrator produces
  finite output with all seven projections.
- The same matrix passes with None, Glitch, Twist, Mobius, and all
  kaleidoscope solids represented in the current lens catalog.
- `After Lens` preserves equality for symmetry-equivalent kaleidoscope inputs.
- Every pair of Planar Warp 1/2 modes, including Projected Vector Noise and
  Projected Curl Flow, is either admitted or rejected by an explicitly listed
  coordinate-semantic rule; no valid pair is collapsed into one stage.
- Every projection and layout/hemisphere suboption matches the normative
  projected-noise matrix independently for Planar Warp 1, Planar Warp 2, and
  Noise Contour (Projected).
- Cut-layout projections do not reject Surface Noise and do not add a
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
- Full-config round trips preserve inactive Lens, Surface Noise, Function, and
  both planar-warp fields while a different dynamic schema is visible, for both
  valid and pending-invalid requested snapshots.
- Every visible enum after each tombstone hole round-trips between dense GUI
  ordinal and stable storage ID; tombstones never appear in a dropdown or new
  export.
- Every rejection path has an exact, actionable warning and no generic
  fallback in expected configurations.

### Compatibility

- Every legacy enum number decodes deterministically.
- Accepted and requested legacy snapshots migrate independently.
- Every two-planar-noise pairing retains both planar stages. Tangent Noise plus
  Projected Vector/Curl retains the projected stages. Tangent Noise plus Stereo
  Noise follows the documented precedence and names the discarded Stereo
  candidate in its import notice.
- Legacy Lens Mix zero and fractional values follow the documented conversion
  and produce a notice.
- New export -> import is lossless for all Surface Noise parameters.
- Old imports produce a visible migration notice and never emit an obsolete
  option on the next export.
- The legacy Noise Contour ID decodes as Noise Contour (Projected); selecting
  Noise Contour (Sphere) exports only its new ID.
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
stack, so that combination is rejected. The retained Noise Contour (Projected)
option supports planar-warp composition under the published projection matrix.
Ignoring the planar warps or silently changing source domain would violate the
no-side-effect contract at the rendered-behavior level.

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
Scale, and bounded coordinate-free vector fields. Retain Projected Vector Noise
and Projected Curl Flow in both planar-warp menus, remove legacy Stereo Noise,
and migrate Tangent Noise out of the Lens menu. Retain Noise Contour (Projected)
and add Noise Contour (Sphere) on the shared scalar field after displacement
migration is stable. Preserve obsolete numeric values only as versioned import
tombstones, with a dense-GUI-to-stable-storage adapter.

This produces more combinations than Stereo Noise while eliminating its
projection dependency. Surface Noise no longer needs projection, pole, seam,
coordinate-bound, or scale-strength rejection rules; explicitly projected
noise retains honest chart-compatibility warnings.
