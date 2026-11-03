# Tangent-field lenses

Status: IMPLEMENTED on 2026-11-02. The design was reviewed against source
77265400d, with owner decisions D1–D6 recorded. Acceptance checks and numerical
clarifications are recorded below.

## Goal

Treat the sphere lenses that tangent fields generate as one family with
shared machinery:

- **Fields without a closed-form flow** share one field interface, one
  closed-form evaluator and one stepping routine. New lenses of this kind
  become parameter sets over that kernel, not new kernels:
  - Rossby–Haurwitz waves;
  - vector spherical harmonic modes;
  - a symmetric flow matching each existing kaleidoscope fold.
- **The Möbius lens** is stored as its generator, a rotation ω and a boost b,
  rather than raw coefficients, so it can be interpolated and never
  degenerates.

Lenses with exact flows keep their own kernels. Twist is unchanged.

Non-goals:

- The symmetry folds (`Kaleidoscope`, `PolyhedralKaleidoscope`,
  `DodecahedralKaleidoscope`, the prisms) and `Glitch`. No field generates
  them: each is non-injective, so no flow produces it. The symmetric flows
  below are designed to pair with the folds, not to replace them.
- Advection: points travelling along a field over elapsed time
  ([Time](#time-and-column-sampling)).
- A general-purpose numerical integrator in the engine.

## Decisions

- **D1 (2026-11-02).** The Möbius lens stores and presents (ω, b),
  everywhere. One vocabulary with the field model, and independent of the
  projection that follows. Cost accepted: pure pans need two sliders moved
  in lockstep.
- **D2 (2026-11-02).** `sphere.lens.mobius.v2` is retired from the engine
  once migrated. The daydream deep-link migrator keeps the v2 names. Old
  MobiusGrid parameter snapshots are rejected by the schema bump.
- **D3 (2026-11-02).** The Rossby–Haurwitz drift ν is authored directly,
  like every other drift clock, and accumulates as a phase.
- **D4 (2026-11-02).** Closed-form Cartesian terms are the only field
  evaluator: general (ℓ, m) modes as fixed polynomials, and one symmetry
  invariant per kaleidoscope fold. There is no runtime Legendre recurrence.
- **D5 (2026-11-02).** Strengths past fold onset are allowed up to the π
  cap. The workbench marks the strength slider beyond a computed fold
  estimate, as a warning and never a clamp.
- **D6 (2026-11-02).** Four terms per field, ℓ ≤ 8. Raising either is an
  operator version bump.

## Baseline at source 77265400d

| Piece | Where | Field view |
|---|---|---|
| `Lens::Twist`, `Lens::AnimatedTwist`, `sphere.lens.twist.v2` | `core/render/pullback/lens.h`, `core/math/lenses.h` | Exact flow of the zonal field θ(p·a)(p × a) |
| `Lens::Mobius`, `sphere.lens.mobius.v2` | `lens.h`, `core/math/mobius.h`, `operators/sphere.h` | Exact flow of a conformal field, stored as raw coefficients a, b, c, d |
| `Surface::CurlNoise`, `Surface::DirectNoise`, `sphere.displace.curl.v2`, `sphere.displace.direct.v2` | `core/render/pullback/surface.h` | One fixed step along a noise tangent field: Euler, midpoint, or two midpoints with `parallel_transport` |
| `Surface::PeriodicRipple`, `sphere.displace.ripple.v2` | `core/render/pullback/surface.h` | Closed-form radial displacement about a centre |
| `math::sphere_exp_map`, `math::sphere_exp_map_half_radian` | `core/math/noise_field.h` | The sphere exponential map; the Taylor form asserts arc ≤ 0.5 rad |
| `math::parallel_transport` | `core/math/spherical.h` | Used by the midpoint step; `HS_CHECK`s within 0.01 rad of antipodal |
| `Animation::MobiusFlow` | `core/animation/params.h` | Already a one-parameter subgroup: a = √s·e^{iθ}, d = e^{−iθ}/√s |
| `Animation::MobiusWarpCircular` | `core/animation/params.h` | Writes b = m·e^{−iθ} into the live params each frame, with a, c and d fixed; θ = easing(progress)·2π, m = scale or a live binding |
| Kaleidoscope folds, `sphere.lens.kaleidoscope.v2` | `lenses.h`, `operators/sphere.h` | Nine symmetries (`KALEIDOSCOPE_SYMMETRY_IDS`) with mirror tables |

Two stage kinds carry these: `Stage::Lens` (direction only) and
`Stage::Displace` (direction plus accumulated `path_length`, which feeds the
`DISPLACEMENT` modulation source). Both are `SphereSample → SphereSample`. A
displacement can sit before or after a lens (`SurfacePlacement`).

Users of the field-type lenses: `GridWorms` (`AnimatedTwist`) and
`MobiusGrid` (`Mobius`, `ANIMATED_MOBIUS = true`). Both are on the Phantasm
playlist (`targets/Phantasm/phantasm_playlist.h`). `MobiusRings`,
`MindSplatter` and `GnomonicStars` use the `Animation::Mobius*` classes on
their own `MobiusParams`, outside the composed pipeline.

### The Möbius coefficient problem

`MobiusLensParams` stores eight floats. Blending two coefficient sets can
cross |ad − bc| = 0, so:

- every coefficient field uses `FieldCurve::SNAP`, and preset interpolation
  holds the old lens and then snaps
  (`interpolate(const MobiusLensParams &, …)`, `composed_parameters.h`);
- runtime writes run a `nondegenerate` check that can refuse a slider value.
  `parameter_write_admitted` reports every LENS refusal as
  `DEGENERATE_WARNING` (`composed_runtime.h`, `operators/sphere.h`);
- column-time sampling follows the animation only when both endpoints were
  generated and share a, c and d (`lens_generated`, `lens_phase_continuous`,
  `composed_interval.h`), and otherwise holds the end lens;
- the animation writes b into the live params, so `ANIMATED_MOBIUS`
  carve-outs keep blends and snaps from fighting it (`composed_runtime.h`).

All four exist because the coefficients are not a vector space and the
animation owns a stored parameter.

## Scope

| Lens | Generator | Treatment |
|---|---|---|
| Killing (rigid rotation) | toroidal ℓ = 1 | No new code. `sphere.rotate.v2` is its flow; inside a field it is an (ℓ = 1) toroidal mode |
| Möbius | toroidal + poloidal ℓ = 1 (sl(2,ℂ), 6 reals) | `Stage::Lens`, exact flow by `mobius_exp` |
| Twist | zonal toroidal, Ψ = ∫θ(h) dh | `Stage::Lens`, unchanged |
| Harmonic modes (general VSH) | up to four toroidal or poloidal (ℓ, m) terms | `Stage::Displace` |
| Rossby–Haurwitz | toroidal (ℓ = 1, m = 0) + (ℓ = R + 1, m = R) | `Stage::Displace`, two mode terms |
| Symmetric flows, one per fold | the fold group's lowest invariant, toroidal (chiral) or poloidal (achiral) | `Stage::Displace`, one invariant term |

The design rule: **a field with a closed-form flow is a lens; a field
without one is a displacement taking one fixed step.** The fixed-step
lenses are what the curl and direct displacements already are. They join
that family, report path length like it, and need no new stage.

## Field model

Every tangent field on the unit sphere splits uniquely (Helmholtz–Hodge;
the sphere has no harmonic 1-forms) as

    u(p) = p × ∇Ψ(p) + Π_p ∇Χ(p),      Π_p v = v − (v·p) p,

with Ψ the stream function (toroidal, divergence-free, area-preserving) and
Χ the potential (poloidal, curl-free, compressive). Ψ and Χ may be extended
off the sphere in any way: two extensions differ by a function vanishing on
the sphere, whose gradient is radial there, and both `p × ·` and `Π_p` drop
radial vectors. So every field below is a polynomial in (x, y, z),
evaluated without trigonometry.

### Conventions

- Axes are the engine's: y is the polar axis, as in `lenses::twist_lens`
  and `mobius_transform`.
- **w = x + i z** and **λ = atan2(z, x)**, the conventions `SHMath`
  (`core/math/spherical_harmonics.h`) and `mobius_transform` already
  share. Every sectoral factor cos^m φ · e^{imλ} is w^m, and the
  stereographic coordinate ζ = (x + i z)/(1 − y) (north pole at ζ = ∞)
  has arg ζ = λ.
- A right-handed rotation by α about +y multiplies x + i z by e^{−iα}: its
  velocity is ŷ × p, and it lowers λ by α. So the Möbius conjugation
  β → β e^{−iθ}, which matches v2's b = m e^{−iθ}, is a right-handed
  rotation of the lens by θ about +y.

Tests pin each sign (T1, T5, T6). Where this section and a test disagree,
the test wins and this section is corrected.

## Field interface

New header `core/math/tangent_field.h`, namespace `math::tangent_field`.
A field is a small value type built in a policy's `prepare`:

```cpp
math::Vector velocity(const math::Vector &p) const;  // tangent at unit p
float speed_bound() const;                           // ≥ max |velocity|
float lipschitz_bound() const;                       // ≥ max ‖∇u‖ (covariant), D5
```

‖∇u‖ is the operator norm of the covariant derivative of u on the sphere,
not of the ambient derivative along the tangent plane; the two differ (for
the tetrahedral toroidal term, 2.484 against 2.530). The `TangentField`
concept checks for these members. Its users are the
stepped displacements, and later the noise displacements (rollout step 4).
No field-level exact-flow concept is defined: Möbius is the only lens that
gains an exact-flow path, and it calls `mobius_exp` directly.

### Conformal generator

Möbius maps act on ζ through SL(2,ℂ). Write the generator as
X = [[α, β], [γ, −α]]. Its flow on ζ is ζ̇ = β + 2αζ − γζ². In sphere terms
it is u = ω × p + Π_p b (rotation ω, boost toward +b̂), with

    α = (b_y − i ω_y) / 2
    β = ((b_x + ω_z) + i (b_z − ω_x)) / 2
    γ = ((b_x − ω_z) − i (b_z + ω_x)) / 2

Rotations come out anti-Hermitian and boosts Hermitian. Both directions of
the conversion, and `velocity` for a generator, live in `tangent_field.h`.
The Möbius lens stores (ω, b) (D1) and converts to (α, β, γ) in `prepare`,
a few adds per frame or column.

**Complex functions.** `math::Complex` (`3dmath.h`) has only + − × ÷ and
conjugate. `csqrt`, `ccosh`, `csinh`, `cacosh` and `clog` are added beside
it as general free functions, since nothing in core provides them, rather
than inside `core/math/mobius.h`.

**Exp** (`math::mobius_exp(X, t)`, new, in `core/math/mobius.h`). X² = δ²I with
δ² = α² + βγ, so

    exp(tX) = cosh(tδ) I + (sinh(tδ)/δ) X

Both cosh and sinh(z)/z are even, so the result does not depend on the
branch of δ. sinh(z)/z switches to its Taylor series below |tδ|² = 1e-4.
The exact exponential has determinant 1. Float coefficient storage loses
absolute determinant precision at strongly conditioned boosts because large
products cancel; tests use a rounding bound scaled by coefficient magnitude
and verify finite unit-sphere maps. Cost: one
complex square root, one `expf`, one `sincosf`, about twenty multiplies.

**Log** (`math::mobius_log(M)`, new, used only for migration and tests).
1. Normalise M by a complex √(ad − bc).
2. Choose the sign with Re tr > 0. On a tie (Re tr = 0, a rotation by
   exactly π), choose the sign with Im tr ≥ 0.
3. Take δ = acosh(tr/2) on the principal branch.
4. Set X = (δ / sinh δ)(M − cosh δ · I), with δ/sinh δ → 1 at δ = 0.

After the sign choice tr/2 ≠ −1, so the only zero of sinh δ is δ = 0.
Elliptic maps use the principal projective branch. At an exactly zero complex
trace, the trace-only tie rule cannot distinguish the two equivalent signs.

## Potential terms

A stepped field holds up to four terms (D6). Each term is a closed-form
scalar T(p) with an analytic gradient, plus a kind (TOROIDAL adds p × ∇T,
POLOIDAL adds Π_p ∇T) and an amplitude A. There is one evaluator. Each term
type supplies its own `gradient(p)`, and the field assembles, sums and
bounds them. There is no runtime recurrence and no special-cased kernel
(D4).

### Normalisation

Every term is scaled so that **amplitude 1 means a peak speed of 1**. The
scale is 1/(1.01 × the measured peak of |∇_S T|), measured on a dense grid.
`speed_bound()` = Σ|A_k| is then a valid bound, tight per term, and the
amplitude means the same thing for every term type. The constants are
generated once by `tools/gen_tangent_field_constants.py` into
`core/math/tangent_field_constants.h`, following the existing generated-LUT
pattern (`tools/gen_gamut_lut.py` → `core/color/gamut_lut.h`, "Generated
by … Do not edit."). They are pinned by T8. For modes, the header stores
only the correction from `SHMath::harmonic_scale` to unit gradient peak;
for invariant terms, the full scale. The margin and T8's window
require the tool to find each peak to within 0.5%, as a relative
tolerance: raw peaks span many decades (mode (8, 8) is about 1.6e7 before
normalisation). The same tool records each unit term's covariant Lipschitz
constant L_k, so `lipschitz_bound()` = Σ|A_k|·L_k (D5).

Unsöld's theorem (Σ_m |∇Y_ℓm|² = ℓ(ℓ+1)(2ℓ+1)/(4π) for orthonormal Y_ℓm)
gives a rigorous but loose bound: measured mode peaks run from 0.71 to 0.33
of it. It stays as a sanity check in T8.

### Harmonic modes

Term `Mode{ℓ ≤ 8, 0 ≤ m ≤ ℓ, φ₀}`:

    T(p) = N · Q_ℓm(y) · Re(c · w^m),      c = e^{−i(φ₀ + mλ₀)}

where λ₀ is the field's accumulated drift phase (D3), N is the term's
scale, and Q_ℓm is `SHMath::reduced_legendre`: P_ℓ^m with its
(1 − y²)^{m/2} factor and P_m^m seed divided out (on the sphere that
factor is |w|^m). `SHMath`'s signed m is the special case φ₀ = 0
(cosine) or φ₀ = π/2 (sine).

Modes are evaluated by `SHMath`, extended in place rather than paralleled:

- **Horner coefficients.** For each (ℓ, m), Q_ℓm is a fixed polynomial in
  y of degree ℓ − m and parity (−1)^{ℓ−m}. A `constexpr` table of its
  coefficients for ℓ ≤ 8 (95 floats, 380 B) is generated at compile time by
  running the existing recurrence on polynomial coefficients, as
  `DodecahedralTable` builds its group, so the recurrence stays the one
  definition. `reduced_legendre` then evaluates by Horner in y² (about
  2(ℓ − m) flops against the recurrence's 9(ℓ − m)), and a new
  `reduced_legendre_derivative` evaluates Q′ from the same table.
- **Phased gradient.** A new `spherical_harmonic_gradient(l, m, p, N, c)`
  returns the R³ gradient, sharing `spherical_harmonic`'s complex-power
  loop:

      ∂x = N Q · Re(c · m w^{m−1})
      ∂z = N Q · Re(c · i m w^{m−1})
      ∂y = N Q′ · Re(c · w^m)

  `spherical_harmonic` stays as the value form with signed m.
- **Scale.** N = `harmonic_scale(l, m)` × the generated gradient-peak
  correction ([Normalisation](#normalisation)), computed in `prepare`.
- **SphericalHarmonics effect.** It (ℓ ≤ 4, on the Phantasm playlist)
  moves onto the Horner path with the rest. Its output is unchanged to
  float rounding, and one device A/B confirms its cost.

The field-level orientation is the quaternion
q = `math::make_rotation`(ŷ, `Vector::from_spherical`(λ_o, φ_o)), applied
as `rotate`(u(`rotate`(p, q⁻¹)), q). It is skipped by a per-frame check
when q is the identity. The azimuthal drift ν is a phase clock shared by
all mode terms.

### Rossby–Haurwitz

The Haurwitz stream function (Williamson et al. 1992, test case 6) on the
unit sphere, in engine coordinates:

    Ψ = −ω sin φ + K cos^R φ sin φ cos(R(λ − λ₀))
      = −ω y + K y Re(e^{−iRλ₀} w^R)

Williamson's longitude increases eastward, along ŷ × p, so it is −λ.
T5 converts accordingly.

It is two toroidal mode terms, (ℓ = 1, m = 0) and (ℓ = R + 1, m = R), with
Q_{R+1,R} ∝ y. R ∈ [1, 7], keeping ℓ ≤ 8. The pattern drifts rigidly:
λ₀ is the accumulated drift phase, advanced by the authored rate ν (D3).

The operator's ω and K are the two terms' amplitudes, with the shared 1.01
normalisation margin. The zonal term's physical rate is ω/1.01.
K does not equal the Williamson K: the unnormalised wave term y·Re(w^R) peaks at
1.000, 1.000, 1.155, 1.299, 1.431, 1.553 and 1.666 for R = 1 … 7, so
K = 1.01 × K_Williamson × that factor. T5's reference applies both factors.

The physical rate ν = (R(3 + R)ω − 2Ω)/((1 + R)(2 + R)), for a planetary
rotation Ω, is documented for reference only. In the lens, ω also adds a
near-rigid azimuthal shift of s·ω to each step.

### Symmetry invariants

One term per kaleidoscope symmetry: the lowest-degree invariant of that
fold's group that is not a function of y alone, written in the fold's own
frame so it aligns with the
existing chamber without an orientation parameter. Each was verified
numerically against the engine's mirror tables (`lenses.h`), with the group
closed from the mirrors:

| Symmetry | Group order | Invariant T | Restriction to the sphere | Term cost |
|---|--:|---|---|--:|
| `azimuthal` | 12 | Re(w⁶), phase 0 | pure ℓ = 6 | mode (6, 6) |
| `tetrahedral` | 24 | ∏ᵢ (aᵢ·p) over the three 2-fold axes aᵢ | pure ℓ = 3 | ≈ 30 flops |
| `octahedral` | 48 | x⁴ + y⁴ + z⁴ | constant + pure ℓ = 4 | ≈ 10 flops |
| `dodecahedral` | 120 | (φ²x² − y²)(φ²y² − z²)(φ²z² − x²) | constant + pure ℓ = 6 | ≈ 30 flops |
| `triangular-prism` … `octagonal-prism` (n = 3, 4, 5, 6, 8) | 4n | Re(wⁿ), phase 0 | pure ℓ = n | mode (n, n) |

- **Tetrahedral axes.** In the engine's tetrahedral frame the 2-fold axes
  are, exactly, (0, −√(2/3), −1/√3), (−1/√2, 1/√6, −1/√3) and
  (−1/√2, −1/√6, 1/√3). They are not written as literals. They are derived
  at compile time from `TETRAHEDRAL_MIRRORS` through the generalised group
  closure ([Reuse](#reuse-and-new-code)): for each half-turn R in the
  closed group, a nonzero column of R + I = 2aaᵀ lies along a. The columns
  are left unnormalised (no constexpr square root is needed), since the
  peak normalisation absorbs the scale. Flipping one axis negates T, and
  −T is equally invariant, so the signs are an orientation convention
  equivalent to A → −A. Gradient by partial products:
  ∇T = Σᵢ aᵢ ∏_{j≠i}(aⱼ·p).
- **Lower-degree invariants that are excluded.** The azimuthal fold also
  preserves y, and the prism folds preserve y². Both are functions of y
  alone, with no azimuthal structure, so they are excluded by the "not a
  function of y alone" rule.
- **Dodecahedral gradient.** With A = φ²x² − y², B = φ²y² − z² and
  C = φ²z² − x²: ∇T = (2x(φ²BC − AB), 2y(φ²AC − BC), 2z(φ²AB − AC)).
  The other cyclic orientation, (x² − φ²y²)…, is not invariant.
- **Where invariance comes from.** Each group has no invariants of lower
  degree, so the restriction is a constant plus one pure harmonic, and
  constants have no gradient.
- **Symmetry by kind.** A poloidal term has the fold's full symmetry,
  u(gp) = g u(p). A toroidal term is chiral, u(gp) = det(g) g u(p), because
  p × ∇ flips under reflection.
- **Implementation.** The azimuthal and prism invariants are ordinary mode
  terms with phase 0, so they add no code. The polyhedral three are the
  only new term types.

**Pairing with a fold.** The image of a displacement D placed before a
fold is source(fold(D(p))). If D is G-equivariant, fold(D(gp)) =
fold(g D(p)) = fold(D(p)), so the kaleidoscope's symmetry survives the
displacement:

- a poloidal term keeps the full symmetry;
- a toroidal term keeps only the rotations, because a reflection reverses
  its swirl, so the fold shows a pinwheel;
- any other field breaks the symmetry.

The symmetric operator therefore takes its symmetry from the same
`KALEIDOSCOPE_SYMMETRY_IDS`, so a chain can share one choice.

## Reuse and new code

No concept here gets a second implementation. Where the engine already has
the piece, the spec uses it or extends it in place.

| Concept | Existing code | Use |
|---|---|---|
| Spherical harmonics (power loop, Legendre, normalisation) | `SHMath`, `core/math/spherical_harmonics.h` | Extended in place: compile-time Horner table, derivative, phased gradient |
| Sphere exponential map, parallel transport | `math::sphere_exp_map`, `sphere_exp_map_half_radian`, `math::parallel_transport` | Reused |
| Euler / midpoint stepping | `Surface::finish_step`, `curl_midpoint_step`, `Integrator` | Lifted onto `TangentField` as `step_field`; CurlNoise moves onto it in step 4 |
| Möbius sphere map | `math::mobius_transform`, `MobiusParams` | Reused as the per-sample kernel and the prepared output |
| Rotation | `math::Quaternion`, `make_rotation`, `rotate`, `Vector::from_spherical` | Reused for field orientation |
| Killing flow | `sphere.rotate.v2` | Reused; no lens |
| Twist flow | `Lens::twist`, `lenses::twist_lens` | Unchanged |
| Reflection-group closure | `lenses::DodecahedralTable::make_group` / `reflect`, hard-wired to `DODECAHEDRAL_MIRRORS` | Generalised to a template over a mirror table; `DodecahedralTable` instantiates it; the tetrahedral axes and T6 use it |
| Fold invariance tests | `expect_chamber_fold`, `tests/test_3dmath.h` | Extended for T6 |
| Kaleidoscope symmetry ids | `KALEIDOSCOPE_SYMMETRY_IDS` | Shared by the symmetric operator |
| Clocks and column-time sampling | `ClockState`, `Spacetime::ClockInterval`, the `SURFACE` phase clock | Reused for the Möbius clock and drift |
| Parameter descriptors | `Field`, `TopologyField`, `TopologyGate` | Reused |
| Generated constants | `tools/gen_gamut_lut.py` pattern | Followed by `tools/gen_tangent_field_constants.py` |

New, because nothing equivalent exists:
- the complex transcendental functions;
- `mobius_exp` and `mobius_log`;
- the (ω, b) ↔ (α, β, γ) conversion;
- the `TangentField` interface;
- the three polyhedral invariant terms;
- the D5 advisory engine-contract entry. A search of the engine and
  daydream found no existing soft-range or slider-mark mechanism.

## Stepping

```cpp
// core/render/pullback/surface.h
template <math::tangent_field::TangentField F, typename IntegratorPolicy>
SurfaceResult step_field(const math::Vector &p, const F &field, float strength,
                         bool path_length_required);
```

This is the existing Euler / Midpoint / Midpoint2 machinery (`finish_step`,
`curl_midpoint_step`) lifted from the noise-specific `curl_field` onto any
`TangentField`. Compiled pipelines pick the integrator at compile time. The
interpreter operators carry it as a runtime enum8, as
`CurlDisplaceParams::integrator` does, and dispatch over the three
instantiations.

**Arc length.**
- If `strength · speed_bound() ≤ 0.5`, use `sphere_exp_map_half_radian`,
  the existing path with its assert.
- Otherwise use `sphere_exp_map`.

The bound and strength are per frame (per column under temporal sampling),
so the branch is uniform across a column. The two maps agree to 1e-7 at arc
0.5, so crossing the threshold does not pop. Nothing clamps the field.
`TangentLimit` saturation adds divergence (documented on
`curl_from_gradient`) and would break the area preservation of toroidal
fields.

**Midpoint transport.** Short steps use `parallel_transport` between their
endpoints. Large steps transport the midpoint velocity back by reversing
the known half-step's great-circle rotation. The initial tangent fixes that
path even at antipodal endpoints, so no joint strength × speed cap is
required. Validation keeps the declared finite scalar and harmonic domains.
Reported path length is traveled arc and may exceed π.

**Fold onset and the advisory (D5).** A single step exp_p(s·u(p)) stops
being one-to-one once s·‖∇u‖ reaches order 1. Past onset, several screen
points sample one source point and the image creases into mirrored pleats,
around vortex centres (toroidal) or sources and sinks (poloidal). At a
crease the source is stretched without limit and aliases, which on the POV
display shows as shimmer along the moving crease.

Review measured the dodecahedral term's onset (unit peak speed, so arc
equals s) at 0.1825 rad (poloidal) and 0.2953 rad (toroidal), well inside
the half-radian domain. With default ranges below onset, the full
exponential map is a safety path that is never reached. For every
single poloidal term tested, onset matched 1/L to three digits. Toroidal
terms fold 1.2–4.6× later for the same L.

**Why 1/L is not a safe advisory.** The step's differential, from Jacobi
fields on the unit sphere, is dF = diag(1, cos ρ) + diag(1, sinc ρ)·s∇u,
in the frame (u/|u|, ⊥), with ρ = s|u|. Curvature costs the perpendicular
row a factor cos ρ, so the flat criterion s‖∇u‖ < 1 is not sufficient. The
tetrahedral poloidal term folds at 0.39787, 1.15% before 1/L = 0.40249.
The sufficient condition for a nonzero Jacobian is s·L < cos(s·B), with
L = `lipschitz_bound()` and B = `speed_bound()`. Since
cos x ≥ 1 − x²/2, it holds whenever s·L + (s·B)²/2 ≤ 1, whose positive root
is closed form:

    s_fold = (√(L² + 2B²) − L) / B²       (→ 1/L as B/L → 0)

Strengths above s_fold are allowed throughout the declared scalar domains. Interpreter operators
report s_fold, computed per frame from the parameters with no iteration.
The value reaches daydream through a new engine-contract entry, and the
workbench draws a warning band on the strength slider beyond it. Nothing
clamps, refuses or changes rendering.

The band is deliberately conservative, not an onset prediction. Sampled
single poloidal terms range from about 1% early for high-degree modes to
27.15% early for degree one. The tetrahedral, octahedral and dodecahedral
poloidal terms are 6.24%, 3.81% and 1.74% early, respectively, with the
normalisation margins included. Toroidal terms and mixes can be substantially
earlier because the derivative bound discards directional cancellation.
The curl displacement has no fold check either; its strength range, ±0.5
with |u| ≤ 1, is set by the half-radian domain.

## Möbius in generator form

### Storage and animation

New family `MobiusGeneratorParams` stores the **authored** generator as a
rotation ω and a boost b: six floats, all `FieldCurve::LERP`, with slider
labels Spin X/Y/Z and Push X/Y/Z (D1). The animation never writes them. The
`LENS` `ClockState` holds the clock (θ, m). The policy's `prepare` composes
the rendered generator by acting on the horizontal parts
H_ω = ω_x + i ω_z and H_b = b_x + i b_z:

    H_ω → m · e^{−iθ} · H_ω,   H_b → m · e^{−iθ} · H_b,   ω_y, b_y unchanged

θ rotates the whole field right-handedly by θ about +y, which in the
generator is the conjugation β → β e^{−iθ}, γ → γ e^{iθ}. m scales the
horizontal parts, which is the generator's off-diagonal part (β, γ). Every
composed generator is valid. `prepare` converts it to (α, β, γ) and calls
`mobius_exp` once per frame, or once per column under temporal sampling.
The per-sample kernel stays `mobius_transform`. An unanimated lens has
clock (0, 1).

What this replaces:

- `start_mobius_animation` becomes `start_lens_phase(scale, duration)`, a
  clock-only driver that publishes (θ, m) with the existing easing and live
  scale binding. It never touches params. Its scale, live-bound or not, is
  the clock magnitude m: a multiplier of the authored horizontal parts,
  spin included, not v2's absolute |b|.
- Deleted:
  - `ANIMATED_MOBIUS` and its carve-outs in `blend_params`, `finish_blend`
    and snaps;
  - `lens_generated` and `ClockState::generated`;
  - `lens_phase_continuous`, whose condition is now always true;
  - `sample_lens_phase` and its `valid() ? sample : lens` fallback;
  - the special-case `interpolate` overload;
  - `nondegenerate`, `DEGENERATE_WARNING` and `COEFFICIENT_LIMIT`;
  - `register_lens_fields` and `MOBIUS_PARAM_NAMES`, since the new family
    registers through its `FIELDS` table.
- `parameter_write_admitted` reports LENS refusals with a range warning in
  place of `DEGENERATE_WARNING`.
- Unchanged: `Animation::MobiusWarp`, `MobiusWarpCircular`,
  `MobiusWarpEvolving` and `MobiusFlow`, which the non-composed effects
  (`MobiusRings`, `MindSplatter`, `GnomonicStars`) keep using on their own
  `MobiusParams`.

**Ranges.** Validation is finiteness and ranges:
- each ω component in [−2π, 2π] (radians of spin);
- each b component in [−6, 6] (b_y = 6 is a zoom of e⁶ linear, e¹² in area,
  at the poles).

Ranges bound the authored values only. The clock's m may push the rendered
horizontal parts beyond them, which is valid because every generator is.

**Authoring consequence of D1.** Pure pans of the stereographic image
(parabolic maps fixing the north pole) lie on the diagonal b_x = ω_z,
b_z = −ω_x with ω_y = b_y = 0, so hand-authoring one moves two sliders in
lockstep.

### MobiusGrid equivalence

Both MobiusGrid presets have c = 0 and identical lens coefficients
(`patterns/mobius_grid.shader.json`). For upper-triangular
M = [[a, b], [0, d]]:

- α = ½ log(a/d) (principal), γ = 0;
- with s = a·e^{−α}, which satisfies s² = ad, the (1,2) entry of exp(X) is
  β sinh α / α = b / s, so β = k·b with the constant k = α / (s sinh α).

For the shipped preset (a = −1.072 + 0.304i, b = 0.416, d = 0.7071):

- **Constants.** α = 0.22739 + 1.43263i and k = 0.01902 − 1.60726i.
- **Authored generator.** β̂ = 0.416·k, γ̂ = 0. Converting through
  B = β̂ + conj(γ̂), W = i(β̂ − conj(γ̂)), b_y = 2 Re α and
  ω_y = −2 Im α gives

      ω = (0.66862, −2.86527, 0.00791),     b = (0.00791, 0.45477, −0.66862)

  That renders the authored snapshot (b = 0.416) at clock (0, 1).
- **Animation.** `MobiusWarpCircular` writes b = scale · e^{−iθ}, so the
  clock-only driver publishes the same θ and m = scale / 0.416. MobiusGrid's
  `start_mobius_animation(1.0f, 160)` becomes
  `start_lens_phase(2.4038f, 160)`.
- **Measured.** The (ω, b) form matches the v2 map to 2.5e-15 (double, 500
  directions, θ over a full turn, three magnitudes). Review measured
  1.06e-6 for float32 `mobius_exp` + `mobius_transform` against the double
  v2 map.

The two presets share their lens, so lerping between them changes nothing.
MobiusGrid's render references are re-recorded for float drift only.

The migration asserts c = 0 for MobiusGrid. The general path serves
workbench documents and deep links: any nondegenerate M through
`mobius_log`, then to (ω, b), at clock (0, 1).

v2 admitted coefficients up to 4 with |det| ≥ 1e-3, so a normalised β or
γ can reach about 126, far outside the v3 ranges. The migrator refuses
such a document with a message naming the lens. It never clamps, because
a clamped lens is a different map.

## Operators, catalog and documents

| Operator | Stage | Params | Notes |
|---|---|---|---|
| `sphere.lens.mobius.v3` | Lens | ω, b (6 × LERP) | Replaces v2 (D2) |
| `sphere.displace.potential.v1` | Displace | strength, integrator, orientation (2), drift ν, 4 × {term type, ℓ, m, kind, A, φ₀} | Term type ∈ {mode, tetrahedral, octahedral, dodecahedral}; ℓ, m and φ₀ are gated live for modes only (`TopologyGate`) |
| `sphere.displace.rossby-haurwitz.v1` | Displace | strength, integrator, ω, K, R, ν | Two mode terms |
| `sphere.displace.symmetric.v1` | Displace | strength, integrator, symmetry (`KALEIDOSCOPE_SYMMETRY_IDS`), kind, amplitude, speed | One invariant term rotating bodily about +y; Speed defaults to zero to preserve fold alignment |

Every displacement's `validate()` refuses m > ℓ, ℓ > 8, R outside [1, 7],
non-finite values, and values outside the declared scalar domains. It never clamps.
Every displacement reports the D5 advisory.

Compiled-pipeline policies: `Lens::MobiusGenerator<Provider>` replaces
`Lens::Mobius<Provider>`, and `Surface::Field<FieldProvider, Integrator>`
serves the three displacements. A composed effect gets a policy only when it
adopts one; interpreter operators come first, for the workbench. The
displacements reuse the `SURFACE` phase clock for drift and per-mode phases.

Touchpoints. Land daydream first, then pin, per the mirror-pin ordering.

- **Engine and catalog:**
  - `operators/sphere.h` and `core/render/pullback/operators/tangent.h`;
  - `scripts/engine_catalog.json` and `tests/data/shader_chain_catalog.json`;
  - `tests/shader_chain/*_cases.h` and `tests/shader_chain/sphere_ops.h`;
  - the 36 Möbius lines in `tests/data/chain_capture_fixtures.jsonl`
    (migrate or re-record).
- **MobiusGrid:**
  - `patterns/mobius_grid.shader.json` and its generated preset header
    (`scripts/generate_composed_presets.mjs`);
  - `MobiusGrid::PARAMETER_SCHEMA_VERSION` bumps from 4 to 5. Snapshot
    restore rejects mismatched versions (`core/control/choreography.h`), so
    old snapshots are rejected (D2). Snapshots are in-memory
    `serialize_parameters` captures with no persistent store, so nothing
    saved is lost.
- **Tests to rewrite or remove:**
  - `tests/composed_effect/temporal_sampling.h`:
    `test_mobius_grid_lens_follows_phase`, and
    `test_mobius_grid_lens_unsupported_transitions` (probably deleted);
  - `case_pullback_mobius_degenerate` in `tests/death/pullback_cases.h`,
    with `DEATH_CASE_COUNT` re-derived from the table.
- **daydream:**
  - the `shader/` mirror;
  - a chain-snapshot migration step, schema 4 → 5, in
    `src/workbench/shader/shader_deeplink.js`, chained after the existing
    3 → 4 step. It rewrites `sphere.lens.mobius.v2` entries and their eight
    coefficient parameters through `mobius_log` into (ω, b);
  - the D5 advisory: a new engine-contract entry and the strength-slider
    warning band, with daydream's engine-contract tests.
- **Indexes:** Holosphere README rows for the new files, including the
  generator tool and the generated terms header.

`sphere.lens.mobius.v2` is removed from the engine once nothing references
it. The daydream migrator is then the only place that still knows the v2
names.

## Time and column sampling

Every lens here is a pure function of (p, parameters at t). Nothing is
integrated over time.

- **Möbius.** L_t(p) = exp(X(t))(p), with X(t) composed from the clock.
  Scaling a generator by t is genuine advection along its field;
  `MobiusFlow` already does this.
- **Stepped fields.** D_t(p) = exp_p(s(t)·u(p; c(t))). Time enters only
  through the drift phases e^{−imλ₀(t)}, closed form in the accumulated
  phase. Orientation and per-term phases φ₀ are static parameters, animated,
  like any parameter, by the existing modulation sources. Under
  `TEMPORAL_SAMPLING`, `prepare` evaluates the phases at each column's time
  through `Spacetime::ClockInterval`, as the twist clocks do. Each sample
  evaluates one field and takes one step. Drift is an accumulated phase,
  not ν·t, so a mid-run change of ν never jumps the pattern.
- **Not provided: advection by a stepped field.** Points do not travel along
  u over time. Φ_t(p) for a field without a closed-form flow needs an ODE
  solve over all elapsed time, which has unbounded cost and drifts
  numerically. A two-phase "streaming" imitation (staggered sawtooth
  strengths blended in colour) is closed form, but costs two source samples
  per pixel and is out of scope.

Symmetric Displace accumulates its Speed in radians per frame, in the range
[-0.1, 0.1]. Its field at phase θ is `R_y(θ) u(R_y(-θ) p)`, rotating the whole
invariant rather than advancing a harmonic phase. This rotates the polyhedral
invariants as well as the sectoral ones. Zero Speed holds the current rotation;
negative Speed reverses it. The existing field phase interval and runtime codec
preserve the rotation across column samples and snapshot restores. Snapshots
from before the Speed control, without that parameter or the instance's runtime
entry, restore at zero rotation.

Temporal continuity: c(t) is continuous, the exponential map is smooth, and
the half-radian switch does not change the value, so D_t is continuous in t
within a frame and across column boundaries.

## Cost

Per-sample flop counts are estimates from the formulas above, not
measurements. "Assembly" is the cross product or projection, the scale and
the half-radian exponential map, about 30 flops. All new kernels are
flash-resident; no ITCM pins are proposed.

| Lens | Term(s) | Per sample, with assembly | Per column (temporal) |
|---|--:|--:|--:|
| Möbius v3 | — | unchanged (`mobius_transform`) | + one `mobius_exp`; measured net preparation overhead ≈ 2,381 cycles per prepared column |
| Mode (ℓ = 6, m = 5) | ≈ 40 | ≈ 70 | one phase factor |
| Rossby–Haurwitz, R = 4 | ≈ 35 | ≈ 65 | one phase factor |
| Symmetric: octahedral | ≈ 10 | ≈ 40 | none |
| Symmetric: tetrahedral or dodecahedral | ≈ 30 | ≈ 60 | none |
| Symmetric: azimuthal or prism n | ≈ 6(n − 1) | ≈ 30 + 6(n − 1) | none |
| Orientation, when not the identity | — | + 2 quaternion rotations (≈ 50) | one quaternion |

At 288 × 144, each 10 flops per sample is roughly 0.7 ms per frame at about
one flop per cycle, against a 62.5 ms frame at 480 RPM. The existing curl
displacement evaluates a 3D simplex gradient per sample, which is expected
to exceed every row above. That expectation is checked when the first
effect adopts a field displacement.

Memory: parameter structs are fixed size (four terms). The flash tables are
`SHMath`'s Horner coefficients (≈ 380 B), the generated peak corrections
and Lipschitz constants (≈ 380 B). No heap, stack or arena tables.

### Device acceptance, 2026-11-02

Matched shipping-profile captures used the same board, compiler, flags and
frame sequence within each comparison:

- SphericalHarmonics, COM4: 218 windows and a full mode cycle. Raster cycles
  per frame fell 4.10% overall; twelve morph groups gained 0.71–1.99% cost,
  at most 0.348 ms/frame. Neither image spilled a frame.
- MobiusGrid, COM4: 68 matched windows covering repeating lens phase cycles.
  Net column preparation added 2,381 cycles/call, or 0.571 ms/frame at the
  observed 144 prepared columns. Total render time rose 0.601 ms/frame
  (1.55%). These are pipeline deltas, not an isolated exponential benchmark;
  the original 300-cycle estimate was low. The capture does not cover a full
  preset cycle and is not a whole-cycle cadence report.
- AshCloud noise consolidation, COM3 and COM4: 1,056 steady frames per board.
  Shader time rose 0.956–0.958% (about 0.44 ms/frame), total render time about
  0.69%, with no steady spills before or after. This is accepted as negligible
  for cadence. Baseline root/wall precision differed by about 10–17 ppm and
  misses the parser's 5-ppm check; matched DWT scope counters carry the cost
  comparison. The shared-path captures pass that precision check.

Phantasm's shipping memory gate measured RAM1 code 196,280 B (+112 B), RAM1
variables 313,760 B (unchanged), and FLASH data 781,228 B (+456 B) against the
baseline. Holosphere and Phantasm both pass size and layout gates.

Windows strict/fast and WASM fast MobiusGrid references are re-recorded.
Linux fast-math exercises the MobiusGrid sequence with a nonblank-output check;
a separate render reference is outside the required validation scope.

## Testing

All native tests join the CTest suite and CI module list.

| ID | Test |
|---|---|
| T1 | For random (ω, b): the central difference (t = ±1e-3) of `mobius_transform(p, mobius_exp(X(ω,b), t))` matches ω × p + Π_p b within 1e-3 on 1,000 directions, including both poles. A negated-ω control must fail |
| T2 | `mobius_log ∘ mobius_exp` is the identity for rotation angle in (−π, π). `mobius_exp ∘ mobius_log` matches M projectively for random, parabolic, identity, near −I, tr = 0 and MobiusGrid inputs |
| T3 | At every corner of the parameter domain, the determinant error obeys a coefficient-scaled float rounding bound and mapped directions remain finite and unit length |
| T4 | Toroidal terms have zero surface divergence and poloidal terms zero surface curl, by finite differences |
| T5 | The Rossby–Haurwitz terms match an independent spherical-coordinate implementation (atan2, cos, double precision) of the Williamson formula, signs included |
| T6 | For each of the nine symmetries: T(`kaleidoscope_lens`(p, symmetry)) = T(p) on random directions. With the group closed from each mirror table, poloidal u(gp) = g u(p) and toroidal u(gp) = det(g) g u(p). A displacement-then-fold image check confirms full symmetry (poloidal) and rotation-only symmetry (toroidal) |
| T7 | `SHMath`'s Horner `reduced_legendre` and its derivative match an independent double-precision recurrence for every (ℓ, m) with ℓ ≤ 8 (libc++ has no `std::assoc_legendre`, so the reference is hand-written). `spherical_harmonic_gradient` matches a finite difference of `spherical_harmonic`, and the existing SphericalHarmonics tests pass unchanged |
| T8 | Every unit-amplitude term's sampled peak speed lies in [0.985, 1.0), so `speed_bound()` is valid and tight. Mode peaks are at or below the Unsöld bound. The generator tool's output matches the committed header |
| T9 | Fold onset per term, for the Rossby–Haurwitz and symmetric configurations, and for random four-term mixes: the smallest s at which the step's Jacobian determinant (Jacobi-field form, checked against finite differences) reaches 0 on the grid. The advisory s_fold lies at or below it in every case. No fixed relative tightness is required; the degree-one poloidal example is 27.15% early. Defaults use strength 0.05 rad, below the sampled single-term onsets |
| T10 | `validate()` refuses m > ℓ, ℓ > 8, R ∉ [1, 7], non-finite values and scalar-domain overages. Every integrator remains finite at the scalar range endpoints; known-path midpoint transport crosses antipodes and reports traveled arcs beyond π. The half-radian assert never fires across the admitted domain |
| T11 | MobiusGrid: v3 and v2 maps agree within 1e-5 on a direction grid for both presets over sampled (θ, m) |
| T12 | Generator interpolation: lerped random presets give finite, continuous maps (a small Δt gives a small map change) |
| T13 | Catalog, shader-document, snapshot-version and daydream deep-link 4 → 5 migration tests, including refusal of a v2 lens whose migrated (ω, b) leaves the v3 ranges |
| T14 | The advisory reaches the engine contract, and the daydream slider shows the band without changing the parameter |

On the tip only: full suite, wasm and wasm-debug smoke, firmware size and
layout gates, the documentation coverage gate, and one re-record of the
MobiusGrid render references.

**Device.** MobiusGrid ships on Phantasm, so one pinned capture after the
Möbius v3 change validates the per-column `mobius_exp` cost. There is no
device work for the displacements until a playlist effect adopts one.

## Rollout

1. **Math, no visible behaviour change.**
   - the complex functions, `mobius_exp`, `mobius_log` and the
     (ω, b) ↔ (α, β, γ) conversions;
   - `SHMath` extended (Horner table, derivative, phased gradient);
   - the group closure generalised;
   - `tangent_field.h` with mode and invariant terms;
   - `tools/gen_tangent_field_constants.py` and its generated header;
   - T1–T9;
   - one device A/B of SphericalHarmonics on the Horner path. Its output is
     unchanged to float rounding and the bar is "negligible or better".
2. **Möbius v3.**
   - generator family, clock-only `start_lens_phase`,
     `Lens::MobiusGenerator` and `sphere.lens.mobius.v3`;
   - the deletions listed above;
   - MobiusGrid migration and schema bump, daydream migrator, catalog;
   - v2 removed; T11–T13 for the lens; device capture.
3. **Field displacements.** `step_field`, `Surface::Field`, the three
   displacement operators in the interpreter and workbench, and the
   advisory; T10, T13 and T14.
4. **Consolidation.** `Surface::CurlNoise` and `DirectNoise` re-expressed
   as `TangentField` types over `step_field`.
   - The curl field type carries CurlNoise's prepared fast path
     (`OWNED_GRADIENT_PREPARE`, `INLINE_UNIT_GRADIENT`) inside its
     `velocity`, so generated code keeps the inlined unit-gradient core.
   - Gated on a two-board device A/B, since CurlNoise is hot in shipping
     effects. The acceptance bar is "negligible".

## Rejected alternatives

- **A generic ODE integrator for every field lens.** N field evaluations per
  sample, iterative, and strictly worse for Twist and Möbius, which have
  exact flows.
- **A runtime Legendre recurrence for modes.** About 4.5× the flops of
  Horner on fixed polynomials, and it still needs an orientation rotation
  to express any fold-aligned invariant (D4). The recurrence survives as
  the compile-time definition of the Horner table.
- **A new harmonic evaluator beside `SHMath`.** A second implementation of
  one concept; `SHMath` is extended instead.
- **Unsöld normalisation.** Rigorous, but it makes "amplitude 1" mean
  anything from 0.33 to 0.71 of peak speed, depending on the mode. It is
  kept as a sanity bound in T8.
- **Storing Möbius as (α, β, γ).** Pans are a single slider there, but the
  vocabulary is tied to one projection (D1).
- **Saturating the field to stay inside the half-radian domain.** A silent
  change of the field that breaks area preservation. The bound-selected full
  exponential map and known-path midpoint transport handle long steps.
- **Capping strength at fold onset.** Onset depends on the whole term mix,
  and any computable bound is conservative, so a cap would also forbid
  smooth settings (D5).
- **A standalone Killing lens.** Identical to `sphere.rotate.v2`.
- **Field-level `ExactFlow` concept and `Sum` combinator.** No user: Möbius
  is the only new exact-flow path, and the one commuting sum (a twist plus a
  rotation about the same axis) is already expressible.
- **Folding the symmetry kaleidoscopes into the field model.** They are not
  field flows; the symmetric terms pair with them instead.

## Open items

None. The implementation uses:

- initial strength 0.05 rad and a published range of ±π, with the full
  amplitude ranges available independently;
- a 1025 × 512 angular grid for constants, checked against independent
  one-dimensional mode peaks to within 0.01% relative;
- `getParameterAdvisories()`, returning current `{name, foldStrength}` entries
  without writing parameters; daydream polls it and shades both strength tails;
- gated enum8 controls for harmonic degree/order and Rossby wave order;
- snapped strength and amplitude endpoints. Orientation, drift, generator components and phases retain
  their specified continuous interpolation.
