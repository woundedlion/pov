# HyperLattice Octet 58 ms — Design Spec

Status: LANDED. Section 3 Tier 1 landed in `f130ee632`; section 4 Tier 2
landed in `1bb522211`. Section 5 fallbacks were not taken. Sections 1-2
record the pre-landing baseline, including historical preset indices and IDs.
The supporting optimization ledger for the preceding arithmetic work is no
longer retained.

Source of truth for the shipped code: `effects/HyperLattice.h`
(shader entry and per-frame preparation), `core/render/sdf/octet_trace.h`
(`trace_3d`, `trace_4d`, `trace_events`), `core/render/sdf/framework.h`
(`OctetEvents`, `OctetEvents4`, `OctetFramework4::edge_query`),
and `core/render/ray/shade.h` (`Appearance`).

## 1. Target and budget

The display interval is 62.5 ms. A frame that renders in more than that waits
for the next flip, so every Octet preset currently displays at 8 fps. The
target is a mean and peak render time of 58 ms or less on the shipping
selective-O3 `profile` image, which leaves the measured 4.9 % ISR share plus a
margin inside one interval and lifts the presets to 16 fps with zero spills.

Each frame evaluates 10,658 samples (the 144×72 quadrant plus its one-pixel
shader margin). Preserve, clear and timeline stepping cost about 0.3 ms, so
the shader has 57.5 ms, which at 600 MHz is **3,236 cycles per sample**.

| Preset | Index | Render mean / peak ms | Cycles per sample | Reduction needed | Report |
|---|---:|---:|---:|---:|---|
| experimental-octet-flight | 2 | 69.841 / 76.566 | 3,819 | 15 % | [preset 3 shipping](../profiles/shipping/profile_hyperlattice_octet_preset3_teensy_2026-09-28.md) |
| experimental-octet-wide-flight | 4 | 73.299 / 76.425 | 4,015 | 20 % | [preset 5 shipping](../profiles/shipping/profile_hyperlattice_octet_preset5_teensy_2026-09-28.md) |
| experimental-octet-4d-slice | 3 | 177.308 / 188.387 | 9,865 | 67 % | [octet shipping](../profiles/shipping/profile_hyperlattice_octet_teensy_2026-09-27.md) |

Two levers are ruled out before starting. The global-O3 twin of preset 3
renders in 68.249 ms against 69.841 shipping, so optimisation level and code
placement are not the gap. The full-roster Phantasm image has about 1 KB of
ITCM left, so nothing moves into ITCM.

## 2. Where the cycles go

### 2.1 Traversal model

The counts below come from a host model of the shipped traversal: the same
four tetrahedral plane families and single-owner pair rule for 3D, the eight
diagonal hyperplane families for 4D, the shipped cell sizes, wire radii, far
distances, near fades and antialiasing footprints, 10,658 random unit
directions per preset, a random lattice-periodic camera position, and the
shipped support test, coverage clamp, opacity and saturation floor. It is a
model, not a device measurement; §6 says how to confirm it on device.

| Preset | Active streams | Candidates per ray | Misses per ray | Contributing per ray | Layers per ray | Max candidates |
|---|---:|---:|---:|---:|---:|---:|
| octet-flight (3D) | 3.00 | 5.78 | 4.86 | 0.92 | 0.92 | 8 |
| octet-wide-flight (3D) | 3.00 | 6.34 | 5.17 | 1.17 | 1.17 | 8 |
| octet-4d-slice | 8.00 | 14.11 | 12.87 | 1.24 | 1.24 | 20 |

A miss is a candidate whose coverage is zero. In 3D the support test rejects
it before the square root; in 4D every candidate runs the full nearest-edge
query, a square root and the footprint divide before its coverage clamps to
zero. No ray in any preset saturates or exhausts the 64-candidate budget, so
the budget caps buy nothing and the early-outs already in place fire only on
the layer path.

### 2.2 What the shipped code does per ray

The supporting disassembly is no longer retained. The recorded inspection of
the pre-landing 3D shader (`shade<false>`, 912 instructions) and the adapter
constructor (483 instructions) gives the following per-ray fixed cost, all of
it spent on values that are constant per frame or true by construction:

- `Interval::valid` and two `Raycast::finite(const Vector&)` calls through
  long-branch veneers on every ray. The direction comes from the trig table
  and is unit by construction; the interval is fixed per frame.
- The `OctetEvents` constructor is an out-of-line call (it carries the hot
  flash attribute, which is `noinline`). Inside it: two `memset` calls for
  value-initialised arrays, twelve `VDIV.F32` (six pair weights, then two per
  stream in `FrameworkPlaneCursor::initialize`, itself another out-of-line
  call), and 48 `VMRS` flag transfers.
- `shade_events` then `memset`s the 68-byte trace result and `memcpy`s the
  72-byte return value on exit.

Per candidate, the generic `trace_events` loop does the following on every
iteration, including the 84 % that are misses:

- Re-tests every stream for finiteness through integer bit checks on values
  that `initialize` and `advance` already guarantee finite.
- Copies a 48-byte `Contribution` (a 64-bit merge identity, a normal and two
  flags the Octet adapter never sets) into the group buffer and validates it.
- Walks the six pairs with an `owner == index && scale != 0` test each, of
  which three run.
- Runs the group-merge bookkeeping even though `GROUP_CAPACITY` is 1.

The 4D shader (`shade<true>`, 1,141 instructions) does the same with eight
streams, and its candidate calls `edge_query<false, true>` out of line: five
`VRINTA`, a two-level sort, parity and tie logic, then `VSQRT`, then the
footprint divide in `contribution()`.

`Appearance::color` branches on the colour mode at runtime, so the DEPTH
shader every preset uses carries the AXIS divide by `feature_count` and its
integer conversion in the hot function.

### 2.3 Why the 4D preset is structurally slower

Every D4 strut direction `(e_i ± e_j)/√2` is orthogonal to four of the eight
hyperplane normals `½(1, ±1, ±1, ±1)`, so every strut lies in a plane of four
different families. `OctetEvents4` has no owner rule: all eight streams are
active, every crossing evaluates the nearest strut of any class to the
crossing point, and the same strut can be composited from several families'
crossings at different points along the ray. The 3D adapter received the
single-owner correction for exactly this; the 4D adapter did not.

## 3. Tier 1: output-preserving changes within oracle tolerance, both domains

These retain output within the existing rendering oracle tolerance; the
rounding change in §3.5 can move a channel by one code value. Modelled
together they remove 30–35 % of the 3D
per-sample cost, which clears the 15 % and 20 % the two 3D presets need with
margin. They are also prerequisites for Tier 2.

### 3.1 A dedicated Octet trace loop

Replace the `trace_events` instantiation with a loop written for this
adapter. The loop keeps the three (3D) or up to four (4D, after §4) active
cursors' `next` and `step` in registers, selects the nearest with two or three
compares, calls the candidate, and composites a hit immediately. The
coincidence rule survives as a single check: a hit whose `t` is within the
relative tolerance of the previous hit keeps the larger coverage and the
first `t`, exactly as the one-slot group did. There is no group buffer, no
`Contribution` copy, no per-candidate validation of fields the adapter never
sets, and no per-iteration finiteness test on stream distances. The three
early terminations stay: interval end, compositor saturation, and the
candidate cap.

This is the change the ledger's rejected trials were reaching for. Grouped
owner lists, the unrolled incidence switch and reciprocal setup were all
local edits inside a 900-instruction function around out-of-line calls, and
their regressions are layout and branch effects. A loop that owns its body is
the unit to measure.

### 3.2 Inline the adapter setup

Make the `OctetEvents` and `OctetEvents4` constructors and
`FrameworkPlaneCursor::initialize` `always_inline` into the shader
instantiation, and construct the cursor and pair arrays without value
initialisation so the two `memset` calls disappear. The shader stays the one
`HS_HOT_FLASH_MEMBER` unit in cached flash; nothing moves into ITCM.

### 3.3 Validate once per frame

`prepare()` already validates the camera, footprint and geometry. Move the
interval check there, and treat the direction as unit by construction from
the trig table. The per-ray veneered calls go away. Fail-fast is preserved at
the boundary where values enter the frame; nothing inside the loop can
produce a non-finite distance once the initialiser has admitted the stream.
`Raycast::finite` and `Interval::valid` remain plain inline helpers; do not
carry optimisation attributes on them, since every `-Os` caller would pay.

### 3.4 No per-ray pair divides

The pair weight is only used in the compare `residual² × scale > support²`.
Store the numerator `a_owner² × spacing²` and the denominator
`a_i² + a_j² + ⅔ a_i a_j` per pair and compare
`residual² × numerator > support² × denominator`. Track the best pair as a
numerator/denominator pair, comparing `n_a × d_b < n_b × d_a`. The one divide
and one square root then run only on contributing candidates, 0.9 per ray.
Cursor setup needs one reciprocal per active stream (`1 / speed` serves both
`next` and `step`) instead of two divides. The ledger's rejected reciprocal
trial replaced divides with the same count of multiplies inside the generic
loop; this removes nine of twelve divides.

### 3.5 Colour mode and premultiplied finish

`Raycast::Appearance` has a single depth colour path, so the shader carries no colour-mode dispatch.
The compositor returns premultiplied colour, so the scan loop does not
un-premultiply with a divide and three rounds and then re-premultiply with
three more. Rounding once instead of twice can move a
channel by one code value; the existing rendering oracle tolerance of two
code values covers it.

## 4. Tier 2 — single-owner struts for D4

This is the only route that brings the 4D preset near budget. It changes the
4D coverage model to the one the 3D adapter already uses: each strut is
evaluated once, by the family that crosses it most directly, and its coverage
comes from the ray-to-line closest approach rather than from the distance at
whichever crossing point happened to be nearest.

### 4.1 Ownership

There are 12 strut classes `u = (e_i + s e_j)/√2`, `i < j`, `s = ±1`. Each
lies in the four families whose normal is orthogonal to it. Per ray, the
owner of a class is the family among those four with the largest
`|n_f · d|`, ties resolved by family index as in 3D. Only owning families keep
active streams; a crossing of family `f` evaluates only the classes `f` owns.
Modelled on the shipped 4D preset:

| | Shipped | With ownership |
|---|---:|---:|
| Active families per ray | 8.00 | 4.00 |
| Candidates per ray | 14.11 | 10.06 |
| Geometry evaluations per ray | 14.11 full edge queries | 35.67 class evaluations |

A class evaluation is about 20 instructions (§4.2); a full edge query is
173 with a sort, five rounds and a square root.

### 4.2 Class evaluation at a crossing

At a crossing of family `f` at parameter `t`, form the normalised point
`q = n + r` with four rounds, as `edge_query` does today. For each owned
class `(i, j, s)`, the nearest strut of that class is the line through the
nearest lattice vertex on that class's coset: for even parity that is `n`,
for odd parity `n` shifted by the unit step that minimises the transverse
distance. The shipped ledger formulas give the squared distances for the best
pair; the per-class closed forms are the same expressions with `(i, j, s)`
fixed rather than selected, and they must be pinned against the exhaustive
24-direction search that the geometry tests already contain.

The offset `v` from the strut line to `q` lies in the hyperplane of `f` and
is orthogonal to `u`, so it lives in a fixed two-dimensional subspace with
orthonormal basis `(w₁, w₂)` per `(f, class)`. Express `v = v₁ w₁ + v₂ w₂`.
The ray-to-line closest approach is then

```text
dist² = |v|² − (v · d)² / (1 − (d · u)²)
v · d = v₁ (w₁ · d) + v₂ (w₂ · d)
```

`(w₁ · d)`, `(w₂ · d)` and `1 − (d · u)²` are per-ray constants of the
owning `(f, class)` pair: 12 triples, two 4D dot products and one
denominator each. Keep the denominator unreduced and compare
`(|v|² × den − (v·d)²) > support² × den`, as in §3.4, so the ray setup has
no divides. The support test rejects a miss before any square root or
footprint divide; those run only on the 1.2 contributing candidates per ray.

### 4.3 What changes visibly

Struts that today are composited from two or more families' crossings are
composited once. Coverage uses closest approach, so a strut passing near the
ray but crossed obliquely gets the same coverage it would from a direct
crossing. Depth, clipping and the footprint still use the owner's crossing
`t`, as in 3D. A ray parallel to a strut has no event for it, as in 3D. The
shipped 4D preview should be re-captured after this lands.

## 5. Fallbacks if the device measurement lands above budget

After Tiers 1 and 2 the 4D preset models at roughly three times faster,
which is at the line, not comfortably under it. The remaining exact levers
are each under 2 %. The honest fallbacks change what the display shows and
are the owner's decision, not the implementer's:

- Drop the one-pixel shader margin for this direct shader if the segmented
  driver's seam contract allows it: 2.8 % of samples.
- Interlace 4D rows across two display windows, or render the 4D preset's
  quadrant at half horizontal resolution and upsample.
- Retune the 4D preset's far distance or wire radius. `far_distance` is the
  multiplier on every candidate count in §2.1.

Do not spend the ITCM reserve, and do not add optimisation attributes to
shared inline helpers; both are measured dead ends recorded in the selective-O3
ledger.

## 6. Measurement protocol

1. Before changing code, capture preset 2 and preset 3 with the deep profile
   enabled (`HS_PROFILE_DEEP=1` in `tools/profile_one.sh`) so the
   `hl_event_step` and `hl_layer_composite` counters confirm the candidates
   and layers per ray in §2.1 on device. Add a miss counter beside them for
   the same run. The deep image is slower; use it for counts, never for
   timing.
2. Land §3 as one A/B against the current shipping baseline, with the
   standard 70-second held-preset capture on both 3D presets and the 45-second
   capture on 4D. Then land §4 as a second A/B. Do not A/B the sub-items of
   §3 individually inside the generic loop; the ledger shows that produces
   noise and regressions.
3. Acceptance: mean and peak render at or below 58 ms on the shipping
   selective-O3 image for both 3D presets, with zero spilled frames and 16 fps
   observed. The 4D preset has the same acceptance; if it misses, record the
   number and choose from §5.
4. Oracles: the 3D changes must pass the existing prepared
   camera rendering comparison at its current tolerances. The 4D ownership
   change alters the coverage model on purpose and needs the same independent
   ray/line cross-product oracle the 3D single-owner correction used, across
   rotations, scale, and the odd-parity vertex shift of every class.
5. Size gates: each change must pass the full-roster
   size/layout gates. The dedicated loop replaces a generic instantiation and
   should shrink the shader; report the FLASH code delta with the timing.

## 7. Reproducing the model

The counts in §2.1 and §4.1 come from a host script that samples 10,658
random unit directions, applies a random rotation and a random
lattice-periodic camera position, and walks the plane crossings with the
shipped spacing (`0.8165 × cell` in 3D, `0.7071 × cell` in 4D), owner rule,
support test `(wire_radius + footprint/2)²`, coverage clamp, opacity
`(1 − t/far)² × cubic(t / near_fade)` and the `10/65535` saturation floor.
Pixel half angle is `π/288`, scaled by the preset's antialiasing strength.
The 4D ownership count assigns each of the 12 classes to the family among
its four with the largest direction dot and counts crossings of owning
families only. Any implementation of §3 or §4 should regenerate these counts
from the device counters rather than from the model.
