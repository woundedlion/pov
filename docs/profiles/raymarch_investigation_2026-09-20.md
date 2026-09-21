# Raymarch: peak-frame investigation and optimization proposals

This is the baseline investigation. See the
[implementation results](raymarch_implementation_2026-09-20.md) for landed
changes, rejected experiments and current peak measurements. Current profile
rankings exclude the pre-publication setup draw.

Snapshot: `e6b82e7457045fb6f8700539801e8c2db506f989`. Teensy 4.0 at 600 MHz,
GCC 15.2.1, real `POVSegmented<288,4,480>` driver and live interrupts.
Captures and ELF/source attestations are archived under
`build/prof/raymarch_investigation_20260920` and `build/prof/artifacts`.
Experiments remain in `C:/work/Holosphere-raymarch-profile` and
`C:/work/Holosphere-raymarch-deep`; production source was not edited.

## Decision

**The current default effect has effectively no recurring frame-time margin.**
The previous README's **62.22 ms** was a peak, leaving just **0.28 ms** before
62.5 ms. This snapshot reaches **62.678 ms** at frame 349 and spills one recurring
frame. I propose **50 ms as the first peak target**, leaving 12.5 ms of headroom.
That is a target, not a measured promise. Average render time is not an acceptance
metric anywhere in this investigation.

Compiler tuning alone does not get there. The strongest next work is to reduce
distance evaluations and the live state of the marcher, then specialize the
canonical torus kernel. Preserve full spatial resolution, the authored AA band,
18-step budget, analytic normals and background-edge compositing.

## Measured peaks

Each standard pass captures 70 s and 1,088 complete frames. All entries below
are maxima, not window means. Baseline and global O3 use COM3; shader-only O3
uses COM4. The deterministic geometry sequence is seeded identically, but
interrupt phasing can vary, so sub-millisecond deltas need repeated A/B captures.

| Image | Recurring peak | Recurring spills | Priming frame |
|---|---:|---:|---:|
| Current selective-O3 | **62.678 ms** | 1/1087 | 108.468 ms |
| Global-O3 reference | 59.044 ms | 0/1087 | 101.827 ms |
| O3 on Raymarch surface/shader region | 60.105 ms | 0/1087 | 104.223 ms |
| Force-inline SDF leaf only | 62.572 ms | 1/1087 | 108.480 ms |
| Force-inline SDF plus transform wrapper | 62.616 ms | 1/1087 | 108.611 ms |
| Reuse the primary distance as the probe seed | 62.225 ms | 0/1087 | 107.669 ms |

The first frame is a different workload: `pov_segmented.h:409` forces full X
coverage to prime both halves before publication. It covers 20,736 pixels;
the recurring quadrant covers 10,368. Its cost belongs to the build/commit
budget. Current standard reports state startup timing separately and exclude
it from runtime statistics and aggregate rankings.

The force-inline experiment is a **rejection**, not a win: approximately 0.06 ms
is within run variation. The full-roster ITCM cost rises 288 B when both levels
are forced. Shader O3 is useful scaffolding, saving 2.573 ms of observed peak
for 528 B additional full-roster ITCM, but it still leaves only 2.395 ms margin.
Its source patch also adds two once-per-frame diagnostic scopes; those compile
away in the roster image. The mathematical rendering algorithm is unchanged,
but framebuffer parity has not been certified for compiler transformations.

## What actually happens in the worst frame

A separate shipping deep capture uses **one-frame windows**, so these are
individual-frame measurements for frame 349. The instrumented render is
64.285 ms; use 62.678 ms from the standard pass for deadline judgments.

| Work in frame 349 | Time | Work count |
|---|---:|---:|
| Primary trace | 29.674 ms | 8,102 rays |
| Occlusion probe | 12.038 ms | 2,728 probes |
| Foreground and background shading, summed | 12.321 ms | 6,520 shades |
| Plotting, summed | 1.910 ms | 6,429 plots |
| Noise table rebuild | 3.114 ms | 3,456 noise samples |
| Timeline | 0.354 ms | once |
| Preserve previous half / clear | 0.224 ms | once |

Remaining time is scan/setup, coverage logic and instrumentation outside these
scopes. `filter_blend` is inside plotting and must not be added again. The
multiple `vol_shade`/`vol_plot` entries have identical labels at different call
sites; this table sums the raw log entries. The generic parser's label-keyed
dictionary retains only the last duplicate, so using it for these totals would
undercount shading.

An independent counter-only capture confirms frame 349's work:

| Counter | Count |
|---|---:|
| All SDF evaluations | **96,494** |
| Cheap-bound returns | 19,514 (20.2%) |
| Positive precise evaluations requiring Lipschitz correction | 76,605 |
| Primary-loop evaluations | 69,393 |
| Probe-loop evaluations | 26,600 |
| Additional refinement evaluations | 501 |
| Primary overrelaxation rewinds | 2,193 |
| Rays reaching the final allowed primary iteration | 494 |

The primary trace takes 8.56 evaluations/ray at this frame; each probe takes
9.75. Rewinds affect 27.1% of primary rays, since a ray disables overrelaxation
after its first rewind. Reaching the final iteration is not proof of a missed
surface, but **494 such rays make a smaller step cap an unsafe quality shortcut**.
In the counter pass another frame reaches 838 final iterations. Frame 349 is
also the highest evaluation-count frame in that pass.

The precise positive path executes three hardware square roots and a reciprocal
divide. Combining the branches gives approximately **269,593 square roots**
at frame 349, before shading. This is a source/codegen-derived operation count,
not an instruction-retirement counter; axis degeneracies and compiler decisions
must be considered for exact dynamic counts.

Interrupts consume roughly 4.8% of the deep peak's elapsed interval. They are
already included in the scope times. DMA packing/submission is far smaller
than tracing; tuning LED submission cannot recover the desired margin.

## Codegen and Cortex-M7 findings

1. **The expensive path already uses hardware FP.** The actual SDF listing has
   `vsqrt.f32`, `vdiv.f32`, and fused multiply/add/subtract; no per-step `sinf`,
   `atan2f`, `powf` or software-double helper. Radial `sqrt(x²+z²)` is shared by
   the bounds test and precise torus distance. Rewriting the source merely to
   share that sqrt would not remove an instruction.
2. **Inlining is contextual.** Full-roster firmware inlines the primary SDF
   path but calls the 360-byte distance function from the occlusion loop and
   refinement. The single-effect profile inlines more of the probe. Forcing
   just `WarpedVolume::distance` inline moves the outlined function into
   `TransformedVolume::distance`; forcing both removes the remaining call but
   does not materially improve the measured peak. Matched compiler and flags
   do not establish identical hot-path machine code.
3. **The pixel body has substantial live state.** The full-roster volume lambda
   saves `d8-d15` and reserves 192 bytes of stack, plus 32 bytes of core-register
   saves: a 288-byte frame. This is register-pressure evidence, not a claim
   that stack capacity itself is exhausted. The scalar trace, probe brackets,
   transforms and three shader paths share this large function.
4. **The per-hit shader crosses optimization and placement boundaries.** The
   baseline full-roster fragment lambda is 508 bytes at flash address
   `0x6003136c`. It saves `d8-d12`, reserves 52 bytes, and calls harmonic,
   normal, reciprocal-root, lighting and palette helpers through several
   flash/ITCM veneers. Two cube-map noise fetches remain out-of-line. Its
   `FunctionRef` thunk adds an indirect-call boundary. These are real calls;
   they are not evidence that every invocation misses the instruction cache.
5. **Memory headroom limits brute-force specialization.** The measured full
   roster has 192,072 B ITCM code, 4,536 B padding to the next 32 KiB bank and
   only 12,800 B for stack/local growth. RAM2 has 4,224 B free. A large ray cache,
   volume texture or another ITCM bank is not a free resource. Profile-only
   images have much more headroom and can conceal this constraint.

The M7 has separate integer/FP resources and supports dual issue; schedule
independent integer indexing around FP dependency chains. It does not provide
NEON-style floating-point ray vectors. Its 32 scalar FP registers make large
software packets liable to spill. Teensy 4.0 also has **hardware double
precision**; treating all doubles as software-emulated would be the wrong
diagnosis. [Arm M7 TRM](https://documentation-service.arm.com/static/5e906b038259fe2368e2a7bb),
[PJRC Teensy 4.0](https://www.pjrc.com/store/teensy40.html).

Do not assume a Quake-style inverse sqrt beats the hardware on this core:
the existing `fast_rsqrt` has two dependent Newton iterations and FP/integer
register transfers. Likewise, do not quote Cortex-M4's divide latency as an
M7 measurement. A replacement needs a target benchmark and an error bound.
No CPI/LSU stall percentages were measured here; instruction dependencies and
spills identify candidates, not quantified stall attribution. The DWT event
counters are only 8 bits, so whole-frame subtraction would silently wrap.

## Proposed work, in priority order

### 1. Reuse the probe seed and simplify ray state

`probe_occluder` starts at `closest_local`, whose distance `closest_d` was just
computed by the primary trace. Pass it through instead of evaluating the same
pure SDF again. Frame 349 offers **2,728 evaluations to remove**, about 2.8% of
its total. This preserves the exact sampled distance, geometry and coverage.
The straightforward prototype reaches **62.225 ms**, a modest 0.453 ms observed
reduction, for 208 B additional full-roster ITCM. This single-pass result needs
repeated A/B confirmation before assigning a precise saving. It is useful
cleanup, not the solution to the headroom problem; loop branching and compiler
decisions clearly limit the benefit.

Then represent progress as scalar ray parameter `t`, retaining `closest_t`
rather than three coordinates. The back-plane test becomes `t > t_exit`;
reconstruct the shading position once. The current loop recomputes a 3D dot
product at every step and stores a three-component best point. The probe
already tracks scalar travel for its parabola while also carrying full point
state. Removing duplicate state offers arithmetic savings **and** room for
better scheduling. Keep the same rewind, hit, graze and step-floor rules;
validate floating-point boundary behavior when changing the update order.

### 2. Compile a canonical torus kernel for the authored twist

All objects are scaled copies with the same proportions and twist. Normalize
the ray once per torus and specialize the common n=2 kernel; dispatch once per
torus to n=0, n=1, n=2 or a generic fallback. Keep the fallback for all live controls.
This removes per-sample twist-loop control, zero-twist tests and reloads of
invariant gate/radius parameters without multiplying the entire renderer by nine.

For n=2, use `sin(2θ)=2xz/(x²+z²)` and `cos(2θ)=(x²-z²)/(x²+z²)`.
The reciprocal can start independently of the radial square root, breaking a
long dependency chain. Transform numerical tolerances and the world-space
step floor consistently when normalizing scale. Check the disassembly: the
value is shorter dependencies and fewer live variables, not simply fewer C++ lines.

Move view and half-light vectors into object space once per torus:
`dot(Q*n, light) = dot(n, inverse(Q)*light)`. This removes a quaternion rotation
from every shade while preserving the lighting model. Keep a compact concrete
shader for this sole volume caller so the compiler can share surface/UV work
and omit alpha sampling whose result Raymarch discards. Apply selective O3 to
that small kernel. Test ITCM and register pressure after every fusion.

### 3. Skip certified empty ray intervals

Start at the ray's actual bounding-sphere entry instead of its front tangent
plane. The screen-plane squared radius is already computed for culling; one
square root supplies entry/exit. Next intersect local cylinder/slab bounds and
subtract a certified empty central-hole interval. This targets the 41.7 ms
primary-plus-probe workload and should help crowded/grazing peak frames.

The more aggressive extension is a tiny hierarchy of conservative capsules
around the sinusoidal torus centerline. Precompute the canonical hierarchy
per twist, reuse it for every scaled/tumbled object, and let it jump over empty
intervals in both primary and background probes. Exact SDF evaluation still
owns hits, normals and final AA.

**The bound must enclose the existing AA distance field, not just the geometric
tube.** Lipschitz-scaled distance can enter the coverage band before raw torus
distance does. Inflate using a proven bound for that metric and preserve
background-graze minima. A geometry-only capsule expanded by `aa_width` is not
automatically safe. Different sample placement under a finite step budget
also needs image validation even when the empty-space bound is correct.

### 4. Replace perpetual noise rebuilding with coherent temporal evaluation

The 24×24×6 table costs 3.114 ms on the critical frame and is rebuilt even though
the default noise phase advances only 0.0002 turns/frame. First precompute the
face-direction normalization factors (one 24×24 grid shared by six faces,
2,304 B as floats) and reuse them. This is a small deterministic optimization.

The aggressive version keeps bracketing temporal tables, updates the next
table in fixed-size slices, and interpolates continuously. Budget a bounded
amount of work **every frame**; rebuilding the whole table every N frames would
make the peak worse. Double buffering adds 3,456 B, so prefer existing scratch
or explicitly budget RAM1. Adjust temporal spacing to speed and scale, and
fall back at rapid live-control changes. Validate hue error after the palette
mapping, where a small noise error can become a larger color error. Temporal
interpolation is a quality-gated proposal, not yet established imperceptible.

### 5. Bigger redesign: analytic mesh interiors with exact SDF boundaries

The surface is directly parameterizable:

```text
x = (R + r cos(v)) cos(u)
y = A sin(nu) + r sin(v)
z = (R + r cos(v)) sin(u)
```

Bake one canonical mesh per twist, transform it per torus, and rasterize
interiors with analytic UVs/normals. Retain exact SDF tracing near silhouettes,
self-occlusion boundaries and ambiguous coverage. Refine geometry by projected
error at the actual LED locations, including polar rows. This could eliminate
most of the 69,393 primary evaluations while concentrating exact work where
the eye notices it. Its saving is unmeasured and depends on interior coverage;
geometry budgets must fit the existing arenas.

This is my preferred bold architectural experiment over lowering resolution
or blindly reusing last-frame pixels. Tight specular motion and touching
torus edges need explicit temporal validation. A tessellation threshold is an
engineering acceptance criterion, not proof of perceptual equivalence.

## Acceptance and next implementation order

Start with probe reuse, scalar ray state, canonical n=2 specialization and
object-space lighting. Then add conservative interval skipping. Keep the
current renderer as an image oracle and the shader-O3 prototype as a small
measured baseline. Do not add the individual timing deltas together; compiler
and cache interactions require a combined measurement.

Acceptance is **worst observed render at or below 50 ms**, zero recurring
spills, and no perceptual regression. Replay identical seeded orientations and
all four segment clips; extend beyond this 70 s capture, and stress all solids,
twists 0–8, fill/AA extremes, crossings, poles, hue seams and fast parameter changes.
Compare linear-color frames, maximum local edge/highlight error and temporal
flicker, then inspect on the actual POV display. Use adversarial geometry
search to generate worst-case replay frames, rather than accepting a low
average or a favorable random seed.

Finally measure the full-roster implementation, since its inlining differs
from the single-effect harness. Require the full firmware memory gate to pass.
The current observations are repeatable evidence about a particular default
trajectory, not a worst-case execution-time proof across the parameter space.

## Validation and artifacts

Both baseline logs and the probe-seed candidate pass `parse_profile.py validate`;
their root cycle totals agree with wall time within 1.3 ppm. Every device run
uses `profile_one.sh`, per-board locking, a unique capture name, source-diff
attestation and archived profile/full-roster ELF images. Baseline and prototype
full-roster `phantasm` builds pass the memory gate. The native suite on the
probe-seed prototype passes 91 applicable tests; one unrelated replay-form test
is skipped (92 registered, zero failures). This is not exhaustive perceptual
validation of the unimplemented architectural proposals.

The standard [shipping](shipping/profile_raymarch_teensy_2026-09-20.md) and
[global-O3](O3/profile_raymarch_teensy_2026-09-20.md) reports and all three
aggregate profile READMEs were refreshed. Previous reports are archived beside
the raw captures. The experimental code has not been landed on the shared branch.
