# GenerativePalette v2 - Design Spec

Status: **IMPLEMENTED.** The recipe path ships in
`core/color/generative_palette.h` (`PaletteRecipe`, `try_compile`,
gamut-relative chroma, hue torsion, `ChromaBasis::PATH_MINIMUM`), the versioned
bridge in `targets/wasm/palette_bindings.h`, and the recipe UI in the daydream
palette tool. The §16 legacy constructors and their `HarmonyType` /
`BrightnessProfile` / `SaturationProfile` enums are removed, so §16 and §17
describe a completed migration, not remaining work.

## 1. Summary

`GenerativePalette` should become a small compiled color path produced from a
larger, setup-time `PaletteRecipe`. The recipe exposes hue, perceptual
lightness, chroma, domain layout, and interpolation as independent controls.
The compiled palette retains three control-space keys so it remains cheap to copy,
sample, bake, and animate with `ColorWipe`.

The default chroma model changes from an absolute, hue-independent OKLCH
ceiling to a fraction of the available sRGB gamut at the sampled lightness and
hue:

```text
h_final = h_path(t) + torsion * (L(t) - 0.5)
C(t) = min(q(t), headroom) * C_max(L(t), h_final)
```

Here `q` is gamut-relative chroma in `[0,1]`, `headroom` defaults to `0.94`,
and `C_max` is the first sRGB gamut boundary on the constant-lightness,
constant-hue OKLCH ray. The existing gamut mapper remains the final safety
net, but ordinary generated palettes should reach it only for floating-point
tolerance, never as their primary saturation mechanism.

This design keeps OKLCH's useful separation of lightness, chroma, and hue,
adds genuinely new palette families such as isolight hue sweeps, constant-
chroma sweeps, tonal monochromes, neutral bridges, and spectral loops, and
avoids an enum explosion of cross-coupled style names.

## 2. Current system and problem statement

The current profile constructor combines four choices:

- `HarmonyType`: triadic, split complementary, complementary, or analogous.
- `BrightnessProfile`: ascending, descending, flat, bell, or cup.
- `SaturationProfile`: pastel, mid, or vibrant.
- `GradientShape`: straight, circular/mirrored, vignette, or falloff.

It resolves those choices into three OKLCH keys. Lightness is mapped into
`[0.12, 0.67]`, and chroma is authored as:

```text
C = (sat / 255) * 0.23 * sin(pi * L)
```

The sine envelope accounts for the broad loss of chroma near black and white,
but it does not account for the much larger hue-dependent variation of the
sRGB gamut. Consequently, many requested colors are outside sRGB. The final
converter correctly reduces chroma while holding OKLCH lightness and hue, but
the mapped color lands on the gamut boundary. One RGB channel then reaches
zero or one, and 8-bit export turns long near-boundary regions into visibly
flat zero/255 runs.

An exhaustive check through the current WASM palette bridge over all 256 base
hues, using straight analogous palettes with flat brightness, found:

| Saturation | Palettes touching an 8-bit bound | Samples touching a bound |
|---|---:|---:|
| Pastel | 0 / 256 | 0 / 65,536 |
| Mid | 113 / 256 | 13,708 / 65,536 |
| Vibrant | 218 / 256 | 40,210 / 65,536 |

This is not a defect in the hue-preserving gamut mapper. It is evidence that
the generator routinely asks the mapper to define the palette.

The current control model has additional limitations:

1. `FLAT` means one fixed high lightness, not a selectable constant-lightness
   palette.
2. Saturation strength and saturation shape are combined into three coarse
   presets.
3. Base hue also acts as the hash seed for profile variation, so changing hue
   can unexpectedly change saturation, lightness, analogous direction, and
   spread.
4. Hue torsion is fixed globally instead of being an explicit aesthetic axis.
5. `CIRCULAR` means a mirrored domain, not a seamless cyclic path.
6. Opposing colors can only travel around cylindrical OKLCH; there is no OKLab
   Cartesian option that deliberately passes through a neutral midpoint.
7. The browser graph compares raw procedural cosine curves with final,
   gamut-mapped and quantized generative output. It cannot show what part of a
   generative path requested unavailable chroma.

## 3. Goals

The redesign must:

1. Make hue, lightness, chroma, and domain layout independently controllable.
2. Make the default generator gamut-aware before RGB conversion.
3. Provide perceptually meaningful controls with stable behavior across hue.
4. Produce substantially more palette families from a compact three-key
   runtime representation.
5. Keep deterministic generation. Changing one explicit control must not
   reshuffle unrelated controls.
6. Preserve exact endpoint behavior, mirror symmetry, loop continuity, and
   smooth `ColorWipe` transitions.
7. Remain suitable for a 600 MHz Teensy 4.0 and the existing arena model.
8. Keep the browser palette tool byte-for-byte aligned with engine output.
9. Preserve legacy palette output until each effect is deliberately migrated.

## 4. Non-goals

This work does not:

- Replace `ProceduralPalette`; cosine RGB palettes remain useful for deliberate
  channel waveforms and controlled over/undershoot.
- Add a general dynamic list of color stops. V2 keeps three unique authored
  keys; loop closure may repeat the first key as a fourth stop.
- Implement display profiles beyond the engine's sRGB/linear-RGB pipeline.
- Promise constant electrical LED power. Constant OKLab lightness, constant
  linear luminance, and constant channel sum are different constraints.
- Make every possible control morphable during `ColorWipe`. Palette topology
  and interpolation policy remain fixed for the subject of a wipe.
- Add four-color harmonies such as square or tetradic before the runtime model
  has four unique keys.

## 5. Terminology

- **OKLab lightness (`L`)**: the perceptual lightness coordinate used by the
  engine. This is the default meaning of "constant brightness" in this spec.
- **Linear luminance (`Y`)**: a weighted sum of linear RGB. A future calibrated
  photometric mode may target it, but it is not the default palette axis.
- **Luma**: a weighted sum of nonlinear RGB. The UI must not use this word for
  OKLab `L`.
- **Absolute chroma (`C`)**: the radius in OKLCH units.
- **Gamut-relative chroma (`q`)**: `C / C_max(L,h)` for `0 < L < 1`. Equal `q`
  means equal use of the available sRGB chroma at each lightness and hue. At
  `L=0` or `L=1`, where the quotient is undefined, `q` remains the authored
  latent intensity control and realized `C` is defined as zero.
- **Headroom**: the maximum allowed `q` for normal generation. Values below one
  keep colors away from a hard RGB boundary.
- **Recipe**: the cold, expressive authoring configuration.
- **Compiled palette**: the small three-key object sampled by `get()` and
  consumed by `BakedPalette`.

## 6. Design principles

### 6.1 Orthogonal axes, named presets as conveniences

Core controls describe color geometry. Names such as `PASTEL`, `VIBRANT`, or
`BELL` are preset builders that assign those controls; they are not separate
runtime algorithms.

### 6.2 Explicit randomness

The compiler is a pure function of a recipe. Random selection belongs in a
separate generator that creates a fully resolved recipe. Its request carries
base hue and variation seed separately, so editing hue does not alter
lightness or chroma.

### 6.3 Generate in gamut, map only as a safeguard

Default modes choose chroma from the gamut boundary before conversion. The
existing hue/lightness-preserving gamut map remains mandatory for custom,
absolute, interpolated, and numerically marginal colors.

### 6.4 Compile richness into three keys

Harmony rules, curves, and seeded variation resolve at construction. The hot
object stores only three control-space keys plus compact evaluation policy. This
preserves the existing allocation-free API and the `ColorWipe` snapshot
budget.

### 6.5 Alternatives considered

**Lower the existing `CHROMA_PEAK`.** This reduces boundary contact but makes
weak-gamut hues dictate saturation indirectly while strong-gamut hues remain
underused. It does not make the saturation control perceptually consistent.

**Keep generating out of gamut and improve the final mapper.** The current
fixed-lightness/fixed-hue chroma reduction is already the right safety behavior
for this engine. A different mapper can change the appearance of boundary
colors, but it cannot recover chroma the display cannot represent, and it
still makes the mapper define much of every vivid palette.

**Use per-channel RGB clipping.** Rejected because it changes hue and
lightness unpredictably and makes animated paths kink when different channels
hit their limits.

**Adopt OKHSV/OKHSL as the whole palette model.** Those spaces are useful for
picker-style controls, but they do not replace the need for independent path,
curve, domain, and morph semantics. Gamut-relative OKLCH supplies the important
available-chroma normalization while retaining the engine's existing
interpolation and gamut infrastructure. An OKHSV/OKHSL input adapter can be
added later without changing the compiled model.

**Make constant absolute chroma the default.** A safe constant must be limited
by the weakest hue/lightness on the path, leaving much of the gamut unused.
It remains available as `PATH_MINIMUM` for users who value exact constancy over
maximum colorfulness.

**Store a baked LUT inside every generated palette.** Rejected because many
palettes coexist in effects such as Dynamo and RingShower. The existing
explicit `BakedPalette` ownership keeps this memory cost visible and
arena-budgeted.

## 7. Public authoring model

The following names are illustrative C++ API, but their separation is a
requirement.

```cpp
enum class HueMode : uint8_t {
  HARMONY,
  SWEEP,
  CUSTOM,
};

enum class PaletteHarmony : uint8_t {
  MONOCHROMATIC,
  ANALOGOUS,
  ACCENTED_ANALOGOUS,
  COMPLEMENTARY,
  SPLIT_COMPLEMENTARY,
  TRIADIC,
};

enum class HueDirection : uint8_t {
  SHORTEST,
  CLOCKWISE,
  COUNTERCLOCKWISE,
};

enum class AxisCurve : uint8_t {
  CONSTANT,
  ASCENDING,
  DESCENDING,
  BELL,
  CUP,
  CUSTOM,
};

enum class ChromaBasis : uint8_t {
  LOCAL_GAMUT,   // q is a fraction of C_max at every sample
  PATH_MINIMUM,  // one literal C safe across the selected path
  ABSOLUTE,      // advanced: keys carry literal OKLCH C
};

enum class ColorPath : uint8_t {
  OKLCH_ARC,
  OKLAB_CARTESIAN,
};

enum class PaletteDomain : uint8_t {
  STRAIGHT,
  MIRROR,
  VIGNETTE,
  FALLOFF,
  LOOP,
};

enum class SegmentEase : uint8_t {
  LINEAR,
  COSINE,
  SMOOTHSTEP,
};
```

The recipe is a POD assembled at setup time:

```cpp
struct HueControls {
  HueMode mode = HueMode::HARMONY;
  PaletteHarmony harmony = PaletteHarmony::ANALOGOUS;
  HueDirection direction = HueDirection::SHORTEST;
  float base_turns = 0.0f;       // [0,1); generation seed is not stored here
  float spread_turns = 0.07f;   // analogous/split spread
  float sweep_turns = 1.0f;     // may exceed one or be negative
  float custom_turns[3] = {};
};

struct AxisControls {
  AxisCurve curve = AxisCurve::CONSTANT;
  float center = 0.62f;
  float range = 0.0f;
  float custom[3] = {};
};

struct ChromaControls {
  AxisCurve curve = AxisCurve::CONSTANT;
  ChromaBasis basis = ChromaBasis::LOCAL_GAMUT;
  float center = 0.62f;          // q for LOCAL_GAMUT/PATH_MINIMUM
  float range = 0.0f;
  float headroom = 0.94f;
  float custom[3] = {};
};

struct PaletteRecipe {
  PaletteDomain domain = PaletteDomain::STRAIGHT;
  SegmentEase easing = SegmentEase::COSINE;
  ColorPath color_path = ColorPath::OKLCH_ARC;
  HueControls hue;
  AxisControls lightness;
  ChromaControls chroma;
  float hue_torsion = 0.0f;      // radians per unit L
  float falloff_start = 0.90f;   // FALLOFF only; strictly between 2/3 and 1
};
```

The default recipe is deliberately moderate: constant `L=0.62`, constant
`q=0.62`, `headroom=0.94`, analogous hue, cosine segment easing, and no hidden
hue torsion. Effects or preset builders opt into stronger structure.

### 7.1 Deterministic validation and canonicalization

Validation runs once before compilation. It processes fields in the table
order below, reports the first error deterministically, and produces a
canonical recipe plus typed adjustment masks. All numeric fields, including
inactive custom fields, must be finite; any NaN or infinity is rejected before
wrapping or clamping. When a grouped rule finds more than one bad field, it
reports the lowest numeric `PaletteRecipeField`; array elements are checked in
index order.

Constants used below are part of the V2 schema:

```text
MAX_SWEEP_TURNS       = 16
MAX_CUSTOM_DELTA      = 16
MAX_CUSTOM_ABS_INPUT  = 4096
MAX_ABSOLUTE_CHROMA   = 1
MAX_ABS_TORSION       = 4*pi radians per unit L
```

| Field | Canonicalization or validation |
|---|---|
| All enum fields | Reject any value not named by its V2 enum. |
| `hue.base_turns` | Wrap with `wrap01`; never clamp. |
| `hue.spread_turns` | Clamp to `[0,0.25]`; negative input becomes zero. |
| `hue.sweep_turns` | Preserve sign; reject `abs(value)>MAX_SWEEP_TURNS`. A loop additionally applies the integer rule in section 8.1. |
| `hue.custom_turns[]` | Reject any input with absolute value above `MAX_CUSTOM_ABS_INPUT`. Subtract `floor(custom[0])` from all three, preserving relative whole turns while putting key A in `[0,1)`. Reject either adjacent absolute delta above `MAX_CUSTOM_DELTA`. |
| `lightness.center`, `lightness.custom[]` | Clamp each to `[0,1]`. |
| `lightness.range` | Clamp to `[0,1]`; negative input becomes zero. |
| Relative-chroma `center`, `custom[]` | For `LOCAL_GAMUT` and `PATH_MINIMUM`, clamp each to `[0,1]`. |
| Relative-chroma `range` | Clamp to `[0,1]`; negative input becomes zero. |
| Absolute-chroma `center`, `custom[]` | Clamp each to `[0,MAX_ABSOLUTE_CHROMA]`; negative chroma becomes zero. |
| Absolute-chroma `range` | Clamp to `[0,MAX_ABSOLUTE_CHROMA]`; negative input becomes zero. |
| `chroma.headroom` | Clamp to `[0,1]` for relative bases. Canonicalize to `1.0` for `ABSOLUTE`, where it is inactive. |
| `hue_torsion` | Preserve sign; reject `abs(value)>MAX_ABS_TORSION`. |
| `falloff_start` | For `FALLOFF`, reject unless strictly inside `(2/3,1)`. For every other domain, canonicalize it to `0.90`. |

After field canonicalization, compilation applies these combination rules:

- Reject `PATH_MINIMUM` with `OKLAB_CARTESIAN`, `VIGNETTE`, or `FALLOFF`.
- Reject a `SWEEP`/`LOOP` whose canonical sweep fails the integer tolerance in
  section 8.1.
- Harmony and custom loops always use the closing rule in section 8.1; there is
  no implicit extra turn.
- A solid or fully achromatic manually authored recipe is valid and is exposed
  as such by diagnostics. The compiler never perturbs user controls to create
  variety.

The embedded API does not throw and does not return a partial palette:

```cpp
enum class PaletteCompileCode : uint8_t {
  OK,
  INVALID_SCHEMA,
  NON_FINITE,
  INVALID_ENUM,
  HUE_LIMIT,
  NON_INTEGER_LOOP_SWEEP,
  INVALID_FALLOFF_START,
  INCOMPATIBLE_OPTIONS,
};

enum class PaletteRecipeField : uint8_t {
  NONE = 0,
  DOMAIN = 1, EASING = 2, COLOR_PATH = 3,
  HUE_MODE = 4, HARMONY = 5, HUE_DIRECTION = 6,
  BASE_TURNS = 7, SPREAD_TURNS = 8, SWEEP_TURNS = 9,
  CUSTOM_TURNS_0 = 10, CUSTOM_TURNS_1 = 11, CUSTOM_TURNS_2 = 12,
  LIGHTNESS_CURVE = 13, LIGHTNESS_CENTER = 14, LIGHTNESS_RANGE = 15,
  LIGHTNESS_CUSTOM_0 = 16, LIGHTNESS_CUSTOM_1 = 17,
  LIGHTNESS_CUSTOM_2 = 18,
  CHROMA_CURVE = 19, CHROMA_BASIS = 20, CHROMA_CENTER = 21,
  CHROMA_RANGE = 22, CHROMA_CUSTOM_0 = 23, CHROMA_CUSTOM_1 = 24,
  CHROMA_CUSTOM_2 = 25, CHROMA_HEADROOM = 26,
  HUE_TORSION = 27, FALLOFF_START = 28, SCHEMA_VERSION = 29,
};

struct PaletteAdjustments {
  uint64_t wrapped_fields;        // bit n corresponds to PaletteRecipeField n
  uint64_t clamped_fields;
  uint64_t canonicalized_fields;
};

struct PaletteCompileStatus {
  PaletteCompileCode code;
  PaletteRecipeField field;       // stable schema field id, or NONE
  PaletteAdjustments adjustments;
};

bool GenerativePalette::try_compile(const PaletteRecipe &input,
                                    GenerativePalette &output,
                                    PaletteRecipe &canonical,
                                    PaletteCompileStatus &status);
```

On failure, `output` and `canonical` are unchanged. On success, `code==OK` and
the three masks identify every field changed by wrapping, clamping, or
inactive-field canonicalization. Internal preset builders may use a checked
`compile()` convenience only after `try_compile()` coverage proves their
recipes valid.

The versioned WASM bridge returns the same status code, stable field id, three
adjustment masks, and canonical recipe. A failure returns no LUT and leaves any
previous module buffer untouched. JavaScript displays the field error; it does
not retry, clamp, or repair independently.

## 8. Control semantics

### 8.1 Hue controls

Hue is stored in turns in authoring APIs and converted to unwrapped radians in
the compiled keys. Unwrapped hue is essential: a full sweep from `h` to
`h + 2pi` must not collapse into two identical wrapped endpoints.

Harmony resolution uses three ordered keys:

| Harmony | Key hues before direction/unwrapping |
|---|---|
| Monochromatic | `h, h, h` |
| Analogous | `h-s, h, h+s` |
| Accented analogous | `h-s, h+s, h+0.5 turn` |
| Complementary | `h, h+0.5 turn, h` |
| Split complementary | `h, h+0.5-s, h+0.5+s` |
| Triadic | `h, h+1/3, h+2/3 turn` |

Positive hue travel is counterclockwise; negative travel is clockwise. Arc
resolution uses these exact turn-space helpers:

```text
wrap01(x) = x - floor(x)                         // [0,1)
r = wrap01(d)

directed_delta(d, SHORTEST):
  r < 0.5  -> r
  r > 0.5  -> r - 1
  r == 0.5 -> (d < 0 ? -0.5 : +0.5)             // signed tie

directed_delta(d, COUNTERCLOCKWISE): r
directed_delta(d, CLOCKWISE):        (r == 0 ? 0 : r - 1)
```

The equality test is exact after the authored inputs have been converted to
the recipe's `float` representation; `0.5` is exactly representable. A zero
wrapped difference remains zero in both forced directions—a forced direction
does not silently manufacture a full turn.

For harmony keys, let `raw[i]` be the base hue plus the signed table offset
above. Compilation resolves:

```text
h[0] = raw[0]
h[i] = h[i-1] + directed_delta(raw[i] - raw[i-1], direction)
```

Consequently the complementary shortest-path tie travels `+0.5` from `h` to
`h+0.5`, then `-0.5` back to `h`; it cannot become an accidental full turn.
Sampling never re-runs arc selection per frame or per segment.

For non-loop shapes, `SWEEP` resolves to:

```text
h0 = base
h1 = base + 0.5 * sweep
h2 = base + sweep
```

This supports partial sweeps, deliberate long arcs, full spectral turns, and
multi-turn helices. Signed `sweep_turns` is already an explicit winding, so its
interaction with `HueDirection` is:

```text
SHORTEST:         sweep = sweep_turns
COUNTERCLOCKWISE: sweep = +abs(sweep_turns)
CLOCKWISE:        sweep = -abs(sweep_turns)
```

For `LOOP`, `n = round(sweep)` must satisfy
`abs(sweep-n) <= SWEEP_INTEGER_EPS`, where `SWEEP_INTEGER_EPS` is exactly
`1e-6`; compilation rejects it otherwise and then stores integer `n`. The keys
resolve to `base`, `base+n/3`, `base+2n/3`, with the repeated endpoint
represented as `base+n`. This distributes the hue travel evenly across all
three loop segments. Here and in the winding formula below, `round` means
nearest integer with exact half ties away from zero.

`base+n` is interpolation and diagnostic state, not an independently converted
endpoint. After clamping the coordinate, a `LOOP` sampler must handle `t==1`
by returning the exact `t==0` `Color4` path before any hue-to-OKLab conversion.
It must not evaluate `cosf(base+2*pi*n)` or `sinf(base+2*pi*n)` for the endpoint.
The unwrapped closing hue is used only for `2/3 <= t < 1`, winding inspection,
and morph compatibility.

`CUSTOM` accepts three already-unwrapped authored turns. `HueDirection` does
not rewrite its first two segment deltas; their signs and whole turns are
authoritative.

For a harmony or custom `LOOP`, the closing key is resolved explicitly rather
than inferred during sampling:

```text
d_close = directed_delta(raw_first - raw_last, direction)  // harmony
d_close = directed_delta(custom[0] - custom[2], direction) // custom
h_close = h[2] + d_close
```

Thus `h_close` is color-equivalent to the first key while the closing direction
is deterministic. Use `SWEEP` when a loop requires an explicit multi-turn
closing winding.

Every resolved segment, including a loop close, stores a measurable winding
class. For resolved delta `d`, define the canonical principal delta with a
counterclockwise half-turn tie:

```text
r = wrap01(d)
principal = (r <= 0.5) ? r : r - 1
winding_class = round(d - principal)             // signed integer
```

The per-segment integer vector is part of morph compatibility. In particular,
`+0.5` has class `0`, `-0.5` has class `-1`, and full turns retain their signed
winding count.

### 8.2 Lightness controls

Lightness is direct OKLab `L`, not HSV value. `center` and `range` resolve to
three nodes:

```text
lo = clamp(center - range/2, 0, 1)
hi = clamp(center + range/2, 0, 1)

CONSTANT:   center, center, center
ASCENDING:  lo,     center, hi
DESCENDING: hi,     center, lo
BELL:       lo,     hi,     lo
CUP:        hi,     lo,     hi
CUSTOM:     custom[0..2]
```

This replaces the fixed `[0.12,0.67]` mapping and makes flat-dark,
flat-midtone, and flat-bright palettes first-class.

Recommended UI range is `[0.05,0.92]`; the API still accepts `[0,1]` so black
and white endpoints remain expressible.

### 8.3 Chroma controls

The same curve vocabulary applies to chroma, but its values depend on the
selected basis.

For both relative bases, `q` is the requested fraction and `headroom` is an
upper limit, both in `[0,1]`:

```text
strength(t) = min(q(t), headroom)
```

Changing the basis therefore does not silently rescale the strength control.

#### LOCAL_GAMUT

The three nodes carry `q` in `[0,1]`. Sampling computes:

```text
C(t) = strength(t) * C_max(L(t), h_final(t))
```

This is the default and the preferred meaning of perceptual color intensity.
It produces approximately consistent gamut use across hue and lightness.

Suggested named strength presets are:

| Preset | q |
|---|---:|
| Muted | 0.35 |
| Balanced | 0.62 |
| Vivid | 0.86 |
| Maximum | 1.00 with explicit headroom 1.00 |

`Maximum` is intentionally allowed to touch RGB boundaries. It must not be
the default.

#### PATH_MINIMUM

This basis produces literal constant or shaped OKLCH chroma without leaving
the gamut. Compilation finds a conservative boundary over the resolved path:

```text
C_path_max = min(C_max(L(t), h_final(t))) over every domain segment
C(t) = strength(t) * C_path_max
```

At the weakest point, `q=0.62` and `headroom=0.94` therefore means 62% of the
available chroma in either relative basis. `PATH_MINIMUM` differs because the
same `C_path_max` is used everywhere, not because `q` is redefined.

The production bound is continuous; a 256-point bake census is validation only
and must not determine `C_path_max`. For each color-bearing OKLCH segment,
including the loop-closing segment, the compiler:

1. Treats eased progress as `u in [0,1]`. `L(u)` and stored path hue are affine
   in `u`, and applying torsion keeps `h_final(u)` affine.
2. Splits the segment at hue wraps and diamond-angle quadrant boundaries.
3. Within each split, enumerates every parameter root at which `L(u)` crosses a
   LUT lightness boundary or `h_final(u)` crosses the inverse diamond-angle ray
   of a LUT angle boundary. Lightness roots are linear. Diamond angle is
   monotone inside a quadrant, so its bucket boundaries are converted to hue
   rays once and the hue roots are linear as well.
4. Sorts those roots, visits the one cell containing each open interval, and
   includes every incident cell at a root or endpoint. Root calculations use
   double precision and outward index nudges so rounding can add a cell but
   cannot omit one. This traverses every closed `(L, angle)` gamut-LUT cell
   intersected by the continuous segment, including all authored stops.
5. Takes the minimum certified chroma floor of those cells.

The current generator's finite subsampling and empirical guard do **not**
certify a cell floor and cannot be used for `PATH_MINIMUM`. Before this mode
ships, the low entry of every LUT cell is regenerated by an offline interval
proof with this contract:

1. The proof domain for cell `(li, ai)` is the exact closed rational box
   `[li/L_steps,(li+1)/L_steps] × [4ai/A_steps,4(ai+1)/A_steps]` in lightness
   and diamond angle.
2. For a candidate quantized chroma `c=k/SCALE`, the prover covers the entire
   radial prism `L × angle × C`, with `C in [0,c]`. Covering the radial interval
   proves the value is before the first gamut exit rather than merely landing
   in a disconnected in-gamut island.
3. The inverse diamond mapping and the exact runtime OKLab-to-linear-RGB
   operation sequence are evaluated with outward-rounded MPFR interval
   arithmetic at a minimum of 128-bit precision. Every source coefficient is
   represented by an interval containing its exact binary32 value, and every
   modeled binary32 operation is widened to include its possible rounding.
4. The prover adaptively bisects lightness, angle, and chroma. A leaf succeeds
   only when all three RGB intervals lie inside the runtime gamut predicate.
   The union of the closed leaves must cover the complete prism. If a leaf
   cannot be proved before the configured depth/width limit, the candidate is
   rejected; uncertainty can only lower the floor.
5. Integer binary search selects the largest proved `k`. Downward integer
   quantization is inherent in `k/SCALE`; cells touching `L=0` or `L=1` store
   zero. The empirical sampled maximum may remain the high bracket entry,
   because underestimating that entry affects refinement quality rather than
   the floor's safety.

The generator emits a proof manifest containing grid dimensions, precision,
subdivision limits, gamut constants, hashes of the transform coefficients and
floor table, and a replayable subdivision certificate. A separate
`--verify-certificate` path re-evaluates every leaf with outward rounding and
fails on any uncovered box or out-of-predicate interval. The certificate is a
release artifact, not firmware data; runtime table layout and lookup cost stay
unchanged. Normal unit tests verify the manifest/table hashes and exercise
dense samples, while the certificate verifier—not finite tests—is the
continuum-wide proof gate.

Authored lightness extrema at exactly zero or one are included explicitly and
force `C_path_max=0` without relying on a LUT cell. This guarantees constant
chroma for arbitrary `get(t)` calls, not only `i/255` samples.

`PATH_MINIMUM` is valid only with `OKLCH_ARC`; a Cartesian OKLab segment does
not have the specified constant-OKLCH-chroma path. It is also invalid for
`VIGNETTE` and `FALLOFF`: interpolation to an injected black stop intentionally
changes chroma and expands the lightness domain beyond the three authored
keys. Those combinations are rejected during recipe validation rather than
claiming a bound that omits part of the sampled domain. This mode is ideal for
scientifically clean constant-lightness/constant-chroma hue sweeps. It will
usually be less saturated than `LOCAL_GAMUT` because the weakest hue limits
the whole path.

#### ABSOLUTE

The nodes carry literal OKLCH `C`. This is required for custom authored keys,
exact color import, and deliberate out-of-gamut authoring. Out-of-gamut samples
use the existing chroma-reduction mapper. Legacy behavior uses its separate
evaluator rather than this basis. The browser labels `ABSOLUTE` as advanced
and shows gamut reduction explicitly.

### 8.4 Hue torsion

Torsion becomes an explicit recipe control:

```text
h_final = h + hue_torsion * (L - 0.5)
```

It is applied before `C_max` is queried. Otherwise chroma would be normalized
against the wrong hue. Default torsion is zero; named presets may use subtle
positive or negative values.

Compiled keys and snapshots always store the **path hue before torsion**. For
an `OKLCH_ARC`, sampling interpolates that path hue and applies torsion exactly
once. Every gamut query uses the resulting final hue:

```text
h_path  = interpolate(key.h_path)
h_final = h_path + hue_torsion * (L - 0.5)
C_max   = gamut_max_chroma(L, h_final)
```

Importing an already-authored OKLCH color into a torsioned palette first
removes torsion, `h_path = h_authored - torsion*(L-0.5)`, so sampling the key
reproduces its authored final hue. Diagnostics expose both path and final hue.

### 8.5 Color path

`OKLCH_ARC` interpolates lightness, chroma control, and already-unwrapped hue
independently. It preserves colorful transitions and is the default.

`OKLAB_CARTESIAN` evaluates the three keys to gamut-safe OKLab and linearly
interpolates `L`, `a`, and `b`. Complementary colors naturally approach gray
instead of rotating around the hue cylinder. Chroma curves define the key
colors, but the Cartesian segment is allowed to reduce chroma between them.
Torsion and local-gamut normalization are applied while resolving each key;
the resulting final-hue OKLab endpoints are then interpolated without applying
torsion again. The UI calls this mode **Neutral bridge**.

#### Achromatic hue resolution

Hue inheritance is defined from control intent, not merely the realized chroma
at one lightness:

```text
LOCAL_GAMUT key is chromatic: q > 0
PATH_MINIMUM/ABSOLUTE key is chromatic: C >= OKLCH_ACHROMATIC_C
injected layout-black stop is chromatic: false
```

A local-gamut key with `q>0` at `L=0` or `L=1` is therefore chromatic for path
resolution even though its instantaneous `C` is zero. For an `OKLCH_ARC`
segment, effective endpoint hues are selected before interpolation:

```text
both endpoints chromatic: interpolate their resolved unwrapped path hues
only left chromatic:       use left hue for the entire segment
only right chromatic:      use right hue for the entire segment
neither chromatic:         use path hue 0
```

Thus a `q=0` key or injected black stop inherits the chromatic neighbor's hue
and cannot seed an arbitrary hue bloom. `OKLAB_CARTESIAN` needs no angular
special case after its keys resolve to `(L,a,b)`; an achromatic endpoint has
`a=b=0`.

Snapshots preserve each authored path-hue float for lossless round-trip, but a
snapshot's `q=0` hue is inert and must never be used as a chromatic endpoint.
For `0 < amount < 1`, a wipe applies the same pair rule to corresponding keys:
both chromatic hues use coherent winding interpolation, one chromatic hue is
held, and two achromatic hues resolve to zero. At exactly `amount=0` or `1`,
the implementation copies the respective snapshot byte-for-byte. This keeps
snapshot endpoint equality while guaranteeing that any newly nonzero chroma
blooms in the chromatic endpoint's hue. Legacy evaluation retains its existing
physical-`C` threshold behavior from section 10.1.

### 8.6 Segment easing

Easing applies to the local segment coordinate after segment selection:

```text
LINEAR:     p
COSINE:     0.5 - 0.5*cos(pi*p)
SMOOTHSTEP: p*p*(3 - 2*p)
```

Cosine is the new recipe default because its zero endpoint slopes remove the
visible derivative kink where three-key gradients meet. Linear remains
available for uniform parameter speed and exact legacy output.

### 8.7 Domain shapes

The shapes specify stop layout and domain remapping, not color generation:

| Shape | Stops | Contract |
|---|---|---|
| Straight | `a@0, b@0.5, c@1` | One pass through all keys |
| Mirror | `a@0, b@0.25, c@0.5`, reflected | Exact symmetry around 0.5 |
| Vignette | `black@0, a@0.1, b@0.5, c@0.9, black@1` | No flat black interval |
| Falloff | `a@0, b@1/3, c@2/3, black@falloff_start, black@1` | Explicit black tail |
| Loop | `a@0, b@1/3, c@2/3, a@1` | Exact color seam; smooth-ease derivative seam |

V2 uses separate `PaletteDomain` and `PaletteHarmony` types so it can add and
name controls without changing the existing ABI. `PaletteDomain` names
`MIRROR` accurately and adds `LOOP`. The legacy `GradientShape` enum
remains unchanged with `STRAIGHT=0`, `CIRCULAR=1`, `VIGNETTE=2`, and
`FALLOFF=3`; existing source and the WASM integer mapping therefore remain
valid. Legacy `HarmonyType` likewise retains its current values and spelling.
The legacy evaluator also retains its current `a,c,b` circular order. New
`PaletteDomain::MIRROR` recipes use the intuitive `a,b,c` order. The browser
may label legacy `CIRCULAR` as "Circular (mirrored)", but the C++ spelling is
not removed during V2 migration.

`falloff_start` must be finite and strictly inside `(2/3,1)`. Compilation
rejects values outside that open interval rather than clamping them into a
zero-width color segment or black tail. The UI range is `[0.70,0.98]`; the API
accepts any representable float in the specified open interval.

For `LOOP`, the repeated final key uses the explicit sweep or harmony/custom
closing rule in section 8.1. A full-turn sweep therefore closes in RGB without
losing its intended hue direction internally. Cosine and smoothstep easing
make both one-sided derivatives zero at the seam; linear easing guarantees
color continuity only.

`BakedPalette::rebake()` recognizes a V2 source whose `loops_domain()` is true.
It samples entries `0..254` normally and copies entry `0`'s RGB and alpha bits
to entry `255`; it does not sample `t=1` independently. Direct sampling and
baking therefore provide bit-identical loop endpoints even when the internal
hue contains many whole turns.

## 9. Generation and randomization

`GenerativePalette::try_compile(recipe,...)` performs no RNG draws. A fully
resolved recipe contains no seed. A separate `PaletteGenerator` creates recipes
from a generation request:

```cpp
struct PaletteGenerationRequest {
  PaletteFamily family;
  bool auto_base_hue;
  float base_turns;
  uint32_t variation_seed;
  uint32_t sequence_index;
};

bool PaletteGenerator::next(const PaletteGenerationRequest &request,
                            PaletteRecipe &canonical_recipe,
                            PaletteCompileStatus &status);
```

Requirements:

1. Reject an unknown `family` or non-finite explicit `base_turns`. When
   `auto_base_hue` is false, wrap `base_turns` with `wrap01` and report that
   adjustment through the shared `wrapped_fields` mask.
2. Automatic base hue uses a seed-derived origin plus the existing full-period
   Weyl step of 157 over the 8-bit hue ring:

   ```text
   hue8 = (hash8(variation_seed, HUE_ORIGIN_SITE) +
           157 * sequence_index) mod 256
   ```

3. Every other randomized field is hashed from
   `(variation_seed, sequence_index, draw_site)` rather than consuming a
   stateful draw stream. Adding a new field must not change old fields.
4. Base hue and `variation_seed` are separate inputs to generation. Moving a
   hue slider on the resolved recipe preserves every non-hue choice.
5. `PaletteFamily` constrains allowed modes and numeric ranges; it does not
   contain hidden per-effect callbacks.
6. Generator families produce canonical valid ranges by construction; they do
   not repair a resolved recipe. A family marked `requires_variety` excludes
   parameter combinations that can produce a solid fill—for example, its
   monochromatic distribution guarantees nonzero lightness or chroma range,
   and its constant-lightness/constant-chroma distribution guarantees hue
   travel. Loop families emit integral sweeps. If a generator bug nevertheless
   produces a recipe rejected by `try_compile`, generation returns that compile
   error and no palette; it does not retry with hidden randomness.
7. Vignette/falloff black regions are explicit. Local-gamut diagnostics mark
   them as layout black rather than low-gamut color; path-minimum recipes reject
   those domains as specified in section 8.3.

Effects may either construct a fixed recipe or request a constrained family.
The browser randomize action resolves and displays the new controls; later
edits operate on that recipe rather than re-running the generator. The engine
must not silently randomize fields the caller supplied explicitly.

## 10. Runtime representation

The compiled palette keeps three **control-space** keys. Their layout remains
three floats per key, but the middle float's meaning is selected by
`ChromaBasis`:

```cpp
struct ControlKey {
  float L;
  union {
    float chroma; // basis-independent internal spelling
    float C;      // compatibility spelling for absolute/legacy callers
    float q;      // LOCAL_GAMUT interpretation
  };
  float h;        // unwrapped path hue before torsion
};

struct Snapshot {
  ControlKey a, b, c;
};

static_assert(sizeof(Snapshot) == 36);
```

The `C` union spelling preserves the existing field layout for legacy source,
but new code must use basis-aware helpers rather than assuming the slot is
always physical OKLCH chroma. `resolved_oklch_key(i)` evaluates one control key
to actual final-hue OKLCH for diagnostics and color import/export.

`LOCAL_GAMUT` stores `q` directly. At `L=0` or `L=1`, realized `C` is zero but
the latent `q` and path hue remain intact in the snapshot. Sampling that exact
key produces black or white; a later `ColorWipe` away from the singular
lightness still has the original intensity available. No `C/C_max` recovery
is performed, so no `0/0` case exists.

Achromatic intent in this basis is determined from stored `q`, not realized
`C`. A key with `q>0` at black or white retains a meaningful latent hue for a
gradient or wipe leaving the endpoint. Only `q=0` makes its hue semantically
irrelevant. Absolute and legacy modes continue to use their physical-chroma
achromatic threshold. Segment sampling and wipes apply the inheritance rule in
section 8.5; they never interpolate from the inert hue of a `q=0` key.

Importing an already-realized black or white color into `LOCAL_GAMUT` cannot
infer a latent intensity and therefore defines `q=0`. A recipe-authored black
or white key retains the recipe's explicit `q`. Black stops injected by
`VIGNETTE` and `FALLOFF` are layout stops rather than canonical keys and do not
erase their neighboring keys' control values.

`PATH_MINIMUM`, `ABSOLUTE`, and the legacy evaluator store literal `C` in the
same slot. Snapshots are only interpreted together with their palette's
evaluation policy; morph compatibility prevents mixing different meanings.

The palette also stores compact policy (`evaluation_model`, `domain`,
`color_path`, `chroma_basis`, `easing`, the per-segment winding-class vector,
`headroom`, `torsion`, `falloff_start`, and the compiled `C_path_max` when
applicable) and derived stop caches. Hue winding records the intentional whole
turns between keys so a spectral or multi-turn path cannot be mistaken for its
wrapped endpoints.

### 10.1 Explicit legacy evaluator

Legacy constructors select a private `LEGACY_SINE` evaluation model rather
than attempting to express legacy behavior through a V2 `ChromaBasis`. That
model retains the existing evaluator unchanged:

```cpp
enum class EvaluationModel : uint8_t {
  V2,
  LEGACY_SINE,
};
```

```text
env_key  = fast_sinf(pi * L_key)
cmax_key = (C_key >= OKLCH_ACHROMATIC_C && env_key > 1e-3)
           ? C_key / env_key
           : 0
cmax(t)  = linear_interpolate(cmax_key)
C(t)     = max(0, cmax(t) * fast_sinf(pi * L(t)))
h_final  = h_path(t) + HUE_TORSION * (L(t) - 0.5)
```

Here `HUE_TORSION` is the current fixed legacy value (`0.35`). Hue interpolation
also retains the current physical-`C` achromatic rules: inherit the chromatic
neighbor's hue when only one stop is gray, and use zero when both are gray. It
also retains the current stop ordering, circular-only cosine easing, gamut
conversion, float helpers, thresholds, and clamp order. Old constructors and
`from_hsv_keys()` set this internal model; V2 recipes cannot request it. The
implementation may keep the legacy evaluator as a separate function until
migration completes. It must not share a V2 branch whose semantics merely
approximate the formula.

The first implementation must keep:

- `sizeof(GenerativePalette::Snapshot) == 36`.
- `sizeof(Animation::ColorWipe) == 104` on wasm32.
- `Animation::ColorWipe <= TimelineEvent::MAX_ANIM_SIZE` on device.
- `sizeof(GenerativePalette) <= 160` bytes unless a measured RAM review raises
  the limit. This matters in particular to Dynamo's 16-entry inline palette
  buffer.

## 11. Gamut boundary API

Refactor the current gamut-bracket machinery to expose:

```cpp
float gamut_max_chroma(float L, float h);
```

It returns the first sRGB exit chroma minus the existing numerical margin. It
must use the same LUT, bracket scan, tolerance, and bisection as
`gamut_clip_preserve_chroma`; there must not be two boundary implementations.
Exact `L=0` and `L=1` return zero before LUT indexing.

The same module exposes the interval-certified per-cell lower-bound read for
`PATH_MINIMUM`. Boundary refinement answers one point; the cell-floor API
answers a closed region. Their numerical margins and gamut predicate are
shared. The offline certificate verifier defined in section 8.3 proves each
closed region; unit tests compare point refinement against the floor as a
regression check but do not claim to prove a continuum by sampling it.

The gamut mapper becomes conceptually:

```text
C_out = min(C_in, gamut_max_chroma(L, h))
```

The generator can then request a fraction of the same boundary. Because the
generated color is already in gamut, conversion requires one OKLab-to-RGB
transform instead of first converting out of gamut, finding the boundary, and
converting again.

`gamut_max_chroma` is the principal performance risk. Most current consumers
already bake a generated palette into 256 entries. Pixel-hot consumers must
continue to use `BakedPalette`; direct use is reserved for setup, sparse
sampling, and bounded per-object lookup.

## 12. ColorWipe behavior

`ColorWipe` continues to snapshot and morph three control-space keys. Palette
evaluation policy belongs to the subject and does not morph.

Two palettes are directly morph-compatible only when these fields match:

- domain shape and stop topology;
- evaluation model, including legacy versus V2;
- color path;
- chroma basis;
- segment easing;
- the complete per-segment winding-class vector, including the loop close;
- headroom, hue torsion, and `falloff_start`, within exact stored equality.

`falloff_start` participates even when only one side would otherwise expose a
different cached stop coordinate. Falloff palettes with unequal values use the
baked output crossfade, which reaches the target topology exactly.

Effects currently create source and target through the same constructor, so
this restriction matches existing use. Add `morph_compatible()` and an
always-on check when a wipe is scheduled.

Before interpolating compatible snapshots, align the target's global hue phase.
The anchor is key A when A is chromatic in both snapshots; otherwise it is the
lowest-index corresponding key chromatic in both. For anchor `j`:

```text
anchor_delta = directed_delta(to.h[j] - from.h[j], SHORTEST)
phase_shift  = from.h[j] + anchor_delta - to.h[j]  // an integer turn
aligned_to.h[i] = to.h[i] + phase_shift

out.h[j] = from.h[j] + amount * anchor_delta
out.h[i] = out.h[j] + lerp(from.h[i] - from.h[j],
                           aligned_to.h[i] - aligned_to.h[j], amount)
```

For every other corresponding key chromatic in both snapshots, the
target-relative offsets—not its arbitrary absolute whole-turn phase—drive the
formula. A key chromatic on only one side uses that side's aligned hue under the
achromatic pair rule. Thus `[0,.25,.5]` and `[10,10.25,10.5]` have zero global
spin. If there is no shared chromatic anchor, all corresponding keys use the
achromatic pair rule from section 8.5. Exact `amount=0` and `1` still copy the
original snapshots byte-for-byte; phase alignment is internal transition
state.

For `LOCAL_GAMUT`, `lerp()` directly interpolates stored `L` and `q`. When both
corresponding keys are chromatic it coherently interpolates their unwrapped
path hue and preserves the compatible winding class rather than making a new
shortest-arc decision. When either key has `q=0`, it uses the achromatic pair
rule in section 8.5. This preserves `q` even when either endpoint has `L=0` or
`L=1`; realized `C` is derived only when the palette is sampled. A wipe between
equally vivid colors therefore does not become dull merely because their
absolute gamut chroma differs.

`PATH_MINIMUM` has an additional compatibility requirement. After global phase
alignment, all three `L` and path-hue keys, the closing hue delta, torsion, and
the compiled `C_path_max` must be bit-equal. Its resolved color geometry is
then identical throughout the wipe and only absolute `C` changes; convex
interpolation of two values bounded by the shared `C_path_max` remains safe.
If any geometry field differs, key morphing is rejected and the transition
uses the crossfade below. V2 does not attempt to certify the two-dimensional
`(path coordinate, wipe amount)` surface in the first implementation.

For compatible `PATH_MINIMUM` and `ABSOLUTE` palettes, `lerp()` interpolates
absolute `C`. `LEGACY_SINE` uses the current absolute-key morph and legacy
evaluator without V2 reinterpretation.

Transitions between incompatible policies use a baked output crossfade, not a
key morph. This distinction is explicit in API names and scheduling.

### 12.1 Incompatible-policy crossfade

`PaletteCrossfade` is a separate animation, not a mode inside `ColorWipe`:

```cpp
struct PaletteTransitionOutput {
  BakedPalette *display_lut;           // exactly one output pointer is non-null
  PaletteCrossfadeView *sparse_view;
  PaletteCoordinateFn coordinate_map;  // identity when null
};

struct PaletteTransitionResult {
  PaletteTransitionKind kind;          // COLOR_WIPE, PALETTE_CROSSFADE, ERROR
  PaletteTransitionError error;        // NONE, OUTPUT_REQUIRED, OUTPUT_ALIASED,
                                       // INVALID_DURATION, SCHEDULE_FAILED
};

PaletteTransitionResult schedule_palette_transition(
    GenerativePalette &subject,
    const GenerativePalette &target,
    PaletteTransitionOutput output,
    int duration,
    EasingFn easing,
    const bool *paused = nullptr);
```

The dispatcher schedules `ColorWipe` only when `morph_compatible()` succeeds.
Otherwise it schedules `PaletteCrossfade`; it never starts a key wipe and then
falls back mid-animation. Invalid output state or duration returns `ERROR` and
schedules nothing.

Ownership and lifetime are:

- `subject` and `target` are caller-owned and must outlive the animation.
  `target` must remain unchanged. `subject` also remains unchanged while the
  crossfade is active, so it is the stable source endpoint.
- Exactly one output pointer is non-null. `display_lut` is the caller's
  already-allocated, non-aliased LUT currently used for rendering `subject`;
  `sparse_view` is a caller-owned non-allocating view used directly by a sparse
  renderer. The selected output and optional coordinate map must outlive the
  animation. The animation borrows them and does not own or allocate endpoint
  LUTs.
- `PaletteCrossfade` owns only references, progress/easing state, and an
  optional coordinate adapter shared by both endpoints. It allocates no arena
  memory and must fit `TimelineEvent::MAX_ANIM_SIZE`.

Add an allocation-free
`BakedPalette::rebake_crossfade(from, to, weight, coordinate_adapter)` operation.
For each existing output entry it samples both palette objects at the same
mapped coordinate and writes the fixed-point RGB/alpha blend directly into the
output arrays. At `weight<=0` it rebakes only `from`; at `weight>=1` it rebakes
only `to`; neither endpoint goes through blend arithmetic. It must not call
the allocating `bake_blend()` API. Callers that sample sparsely instead of
owning a LUT use a non-owning `PaletteCrossfadeView` over the same two endpoints
and current weight. The animation updates the caller-owned view in place.

On the completion step, the animation first renders the exact target endpoint,
then copies the target's complete compiled state—control keys, evaluation
model, domain/topology, chroma policy, winding vector, torsion, falloff, and
derived caches—into the existing `subject` object. The subject's address does
not change. Subsequent rebakes and direct samples therefore use the target
policy without depending on the animation or target reference. A sparse view
is rebound to the adopted subject on completion so it also releases the target
reference.

Cancellation before completion restores `display_lut` from the unchanged
subject, or rebinds the sparse view to that subject, and leaves the subject
policy untouched. Pausing changes neither the weight nor output. Supplying zero
or two output pointers returns `OUTPUT_REQUIRED`; an aliased LUT returns
`OUTPUT_ALIASED`. Neither error mutates the subject or output, and neither path
allocates hidden storage. This contract adds no persistent arena allocation.

## 13. Reference palette families

These are recipe presets for the tool and examples for effects, not new core
enums:

| Family | Hue | Lightness | Chroma | Path/shape |
|---|---|---|---|---|
| Balanced analogous | Analogous | Constant 0.62 | Local gamut 0.62 | OKLCH/straight |
| Isolight harmony | Triadic or split | Constant selectable L | Local gamut 0.72 | OKLCH/straight |
| Constant-chroma sweep | Full/partial sweep | Constant L | Path minimum | OKLCH/straight |
| Tonal monochrome | Monochromatic | Asc/desc/bell | Bell local gamut | OKLCH/straight |
| Chroma breath | Fixed/analogous | Constant L | Bell or cup | OKLCH/mirror |
| Neutral complement | Complementary | Constant or cup | Local gamut | OKLab/straight |
| Spectral loop | Full-turn sweep | Constant or subtle bell | Local gamut | OKLCH/loop |
| Helix | Full/multi-turn sweep | Ascending/descending | Constant local gamut | OKLCH/loop or straight |
| Vignette glow | Analogous/accented | Bell | Bell local gamut | OKLCH/vignette |

The preset list demonstrates independent axes. Users can start from a family
and then edit any control without the preset being a separate locked mode.

## 14. Browser tool redesign

The browser implementation is the Daydream page `tools/palettes.html` and its
existing companion modules `palette_controls.js`, `palette_math.js`, and
`palette_canvas.js`. Reorganize its generative tab into the following sections.
The procedural tab and its named-palette gallery remain behaviorally unchanged.

### 14.1 Hue

- Base hue.
- Mode: harmony, sweep, custom.
- Harmony when applicable.
- Spread or sweep turns.
- Direction.
- Variation seed and a separate randomize button. Randomization populates the
  visible controls; changing base hue afterward does not randomize again.

### 14.2 Lightness

- Curve.
- Center.
- Range.
- Three key controls when custom.

The label is **OKLab Lightness**, with a tooltip distinguishing it from RGB
luma and LED power.

### 14.3 Chroma

- Basis: local gamut, constant/path minimum, absolute.
- Curve.
- Intensity/center.
- Range.
- Headroom under an advanced disclosure.

### 14.4 Path

- Domain shape.
- Color path: colorful OKLCH arc or neutral OKLab bridge.
- Segment easing.
- Hue torsion.
- Falloff start when relevant.

### 14.5 Diagnostics

The current "RGB Wave Functions" title becomes **Palette Curves** on the
generative tab. The graph supports overlays:

- Final sRGB output, after the same quantization the device preview uses.
- `L(t)`.
- Absolute `C(t)` and boundary `C_max(t)`.
- Gamut fraction `q(t)`.
- Markers where the fallback gamut mapper changed the requested color.

Procedural mode continues to show raw unclamped cosine curves. Generative mode
must not imply that its final byte LUT is an unclamped source function.

The color strip remains the authoritative device preview.

### 14.6 Page structure and legacy coexistence

Change the page title from **Procedural Color Palette Generator** to **Palette
Designer**. Keep the top-level Procedural and Generative tabs. The Generative
tab contains a model switch:

- **V2 Recipe** is the default and exposes the controls in sections 14.1-14.5.
- **Legacy Profiles** preserves the current Base Hue, Harmony, Saturation,
  Brightness, and Gradient Shape controls, their DOM ids, legacy WASM
  `bakeLut` call, and legacy C++ export for regression comparison.

V2 controls use semantic `<fieldset>`/`<legend>` groups rather than a fixed
five-column row. The layout is one column on narrow screens, two on medium, and
three on wide screens. Conditional controls use the `hidden` attribute and
update `aria-expanded`/`aria-controls`; merely graying a control is not enough
because inactive values must not remain in the keyboard order.

Add a V2 preset gallery for the reference families in section 13. Loading a
preset writes every visible recipe field, clears any transition snapshots, and
then compiles normally. It does not create a hidden locked preset mode. The
existing named procedural gallery remains procedural-only.

### 14.7 State and update pipeline

V2 state is a plain object owned by `palette_controls.js`, not reconstructed by
ad hoc DOM reads in `updatePalette()`:

```text
GenerativeToolState {
  schema_version,
  model,                         // V2_RECIPE or LEGACY_PROFILES
  draft_recipe,
  canonical_recipe?,
  generation_request,
  adjustments,
  validation_error?,
  overlays,
  last_good_lut?,
  last_good_diagnostics?,
  transition_from?,
  transition_to?,
}
```

Control events dispatch typed state actions. One render function projects state
back to DOM values, visibility, accessible descriptions, preview canvases, and
export text. Programmatic application of a canonical recipe is guarded so it
does not recursively dispatch input events.

Keep the existing animation-frame coalescing: any number of edits in one frame
causes at most one WASM compile/inspect and one canvas repaint. Each edit has a
monotonic revision; a result is applied only to the revision that requested it.
This makes the contract safe if the bridge later moves to a worker even though
the initial WASM call is synchronous.

`palette_math.js` must not implement V2 harmony, curve, gamut, validation,
transition, or randomization math. It retains the procedural and legacy-profile
implementations, and adds only V2 wire/result adapters, copied-LUT sampling, and
C++/JSON serialization of an engine-returned canonical recipe. Every typed
array view returned from WASM is copied before another bridge call because it
aliases module-owned buffers.

### 14.8 Validation and control behavior

Every V2 edit sends the complete draft recipe to `inspectV2`. On success:

- apply the engine's canonical recipe to state and controls;
- display non-error adjustment chips from the corresponding wrapped, clamped,
  or canonicalized field masks;
- replace the last-good LUT and diagnostics atomically; and
- enable C++/JSON export and transition capture.

On failure, retain the last-good preview but mark it **Last valid preview**,
show the shared error code and field-specific message in an `aria-live` status,
set `aria-invalid` and `aria-describedby` on the responsible control, and
disable export and transition capture. A button in the error status focuses the
field. The page never repairs the recipe in JavaScript and never exports the
last-good recipe as if it matched invalid controls.

The UI proactively hides or disables impossible choices—for example,
`PATH_MINIMUM` with Cartesian, vignette, or falloff—but still submits the full
recipe and treats the engine as authoritative. Numeric controls pair a range
input with an editable numeric input so wrapping, clamping, rejection, and the
limits in section 7.1 are inspectable rather than being silently constrained by
HTML slider bounds.

Control-specific requirements are:

- Harmony appears only for `HARMONY`; spread appears only for harmonies that use
  it. Sweep and direction appear for `SWEEP`; three turn inputs appear for
  `CUSTOM`.
- Chroma labels its value `q` for relative bases and `C` for `ABSOLUTE`.
  Headroom remains visible under Advanced for both relative bases and is hidden
  for absolute chroma.
- `falloff_start` appears only for `FALLOFF`. Loop sweep validation displays the
  exact integer-tolerance error from the engine.
- The existing strip drag-to-zoom continues to rewrite procedural cosine
  parameters only. In V2 it changes the displayed coordinate window without
  changing the recipe; Reset Zoom restores `[0,1]`.

### 14.9 Diagnostics rendering

Replace the generative use of `drawWaveGraph()` with
`drawPaletteCurves(diagnostics, overlays)`. Procedural mode continues to use the
current RGB wave painter. The V2 graph has a legend and two labeled axes:

- the left `[0,1]` axis carries final sRGB, `L`, and `q`;
- the right chroma axis carries `C` and `C_max`, scaled to the maximum visible
  `C_max` and labeled with its numeric range.

Fallback-map samples are vertical markers, not another interpolated curve.
Authored key/stop positions are optional markers. Hover or keyboard movement
over the graph reports `t`, RGB bytes, `L`, `C`, `q`, `C_max`, `h_path`,
`h_final`, and fallback status from the nearest engine diagnostic sample.

For `LOOP`, show a **Bit-exact seam** badge only when LUT entry 255 equals entry
0 in all channels. The strip painter always consumes the copied engine LUT and
must not rebuild entry 255 from CSS interpolation. Diagnostics may show the
unwrapped closing `h_path`; the output color at `t=1` remains the aliased first
entry.

### 14.10 Generation and transition preview

The Randomize button calls the engine's `generateV2` with the visible
`variation_seed`, `sequence_index`, family, and base-hue policy. It replaces the
draft with the returned canonical recipe. Editing base hue afterward changes
only hue; Randomize is the only action that advances `sequence_index`.

Add a **Transition Preview** disclosure with **Capture From**, **Capture To**,
duration, easing, Play/Pause, and Reset. Captures are deep copies of canonical
recipes, not references to the live draft. Each preview frame calls the V2
transition bridge; JavaScript does not reproduce phase alignment,
`PATH_MINIMUM` eligibility, or crossfade dispatch. Display the engine-returned
kind (**Key wipe** or **Output crossfade**) and incompatibility reason. At
progress zero and one, the preview LUT must byte-match the captured endpoint
LUTs. Editing the live recipe does not mutate either capture.

### 14.11 Export and file boundaries

V2 C++ export serializes the complete engine-returned canonical recipe,
including schema-relevant defaults, and demonstrates checked
`try_compile(...)` handling. It never emits a legacy constructor. Add a JSON
recipe export containing `schema_version` and the canonical fields; importing
unknown versions or enum values reports the shared validation error and changes
no current state. Legacy mode retains its existing constructor export.

Implementation responsibilities are:

| File | V2 responsibility |
|---|---|
| `tools/palettes.html` | Semantic markup, responsive groups, thin event wiring, status/transition controls, and accessibility. |
| `tools/palette_controls.js` | Default state, reducer/actions, control metadata, conditional visibility, canonical-result application, and transition captures. |
| `tools/palette_math.js` | WASM adapters, immediate buffer copies, LUT sampling, and canonical C++/JSON serialization; no V2 color algorithm. |
| `tools/palette_canvas.js` | LUT strip rendering and the dual-axis diagnostic graph, markers, legend, hit testing, and seam status inputs. |
| `tests/palette_controls.test.js` | Reducer determinism, visibility, adjustment/error mapping, revision rejection, and immutable captures. |
| `tests/palette_math.test.js` | Buffer-alias copying, result decoding, LUT endpoints, and exact canonical exports using a fake bridge. |
| `tests/palette_canvas.test.js` | Overlay scales, fallback/key markers, hit testing, and loop seam rendering. |

Add one browser smoke test that switches among Procedural, V2 Recipe, and
Legacy Profiles; exercises keyboard access and a validation error; verifies a
successful compile/export; and previews both a key wipe and a forced crossfade.

## 15. WASM contract

The C++ engine remains the source of truth for recipe compilation and palette
sampling. JavaScript owns UI state, not a second implementation of the color
algorithm.

Add a versioned recipe bridge, for example:

```text
PaletteOps.compileAndBakeV2(recipe) -> CompileResult { status, canonicalRecipe, lut? }
PaletteOps.inspectV2(recipe) -> InspectResult { status, canonicalRecipe, lut?, diagnostics? }
PaletteOps.generateV2(request) -> GenerateResult { status, canonicalRecipe? }
PaletteOps.bakeTransitionFrameV2(from, to, progress)
    -> TransitionResult { status, kind, incompatibilityReason, lut? }
```

The wire representation may be a validated embind value object or a packed
versioned POD. It must include a schema version and reject unknown enum values.
Status codes, field ids, and all three adjustment masks exactly match section
7.1. A failed call contains no LUT/diagnostics and does not modify the previous
module buffers. A successful returned LUT keeps the current module-lifetime
aliasing contract.

`bakeTransitionFrameV2` is stateless tool support over canonical endpoint
recipes. It runs the same compatibility, global phase alignment,
`PATH_MINIMUM`, endpoint alias, and output-crossfade rules as the runtime
animations. Its progress is clamped to `[0,1]`; endpoint results are copied from
the compiled endpoint LUTs. `generateV2` uses the exact generator and hash
inputs from section 9. Neither operation delegates color decisions to
JavaScript.

The existing `PaletteOps.bakeLut` entry point and its `GradientShape` integer
ABI remain unchanged: `STRAIGHT=0`, `CIRCULAR=1`, `VIGNETTE=2`, and
`FALLOFF=3`. V2's `PaletteDomain` is carried only by the new versioned bridge;
`LOOP=4` is never passed through the legacy cast. Existing harmony,
brightness-profile, and saturation-profile integer mappings used by legacy
WASM exports are unchanged as well.

`inspectV2` may reuse fixed module buffers and is tool-only. It exposes the
engine's resolved `L`, `C`, `h_path`, `h_final`, `q`, `C_max`, and fallback-map
flag for 256 samples so the graph cannot drift from device math.

## 16. Legacy compatibility and migration

Existing constructors remain available during migration:

```cpp
GenerativePalette();
GenerativePalette(GradientShape, const CPixel &, const CPixel &, const CPixel &);
GenerativePalette(GradientShape, const Pixel &, const Pixel &, const Pixel &);
GenerativePalette(GradientShape, HarmonyType,
                  BrightnessProfile, SaturationProfile, int seed);
GenerativePalette::from_hsv_keys(...);
```

They select the explicit private `LEGACY_SINE` evaluator described in section
10.1, which preserves:

- current HSV-like key resolution and random draw sites;
- current `CHROMA_PEAK`, sine envelope, lightness band, and fixed torsion;
- current linear easing except for circular smoothing;
- current circular `a,c,b` ordering;
- current byte-level output and golden tests.

New code uses:

```cpp
GenerativePalette::try_compile(const PaletteRecipe &recipe, ...);
```

Do not silently redirect old constructors to new defaults. Each effect must be
visually reviewed and migrated intentionally. After all production call sites
and the browser exporter use V2, legacy APIs may be removed in a separate
change.

Generated C++ snippets from the browser include the full V2 recipe or a named
recipe builder. They must never export a constructor whose defaults are needed
to reproduce the preview.

## 17. Implementation phases

### Phase 0 - Baseline

- Commit the boundary-incidence census used in this spec as a repeatable host
  or WASM test utility.
- Record current LUT goldens, `sizeof` values, device cycles for a 256-entry
  bake, and ELF/ITCM size.

### Phase 1 - Shared gamut boundary

- Extract `gamut_max_chroma(L,h)` from the existing mapper.
- Prove the mapper is byte-identical before and after the refactor.
- Add the outward-rounded interval prover and independent certificate replay.
- Regenerate the LUT low entries and commit their proof manifest; the current
  sampled-and-guarded floors are not eligible for `PATH_MINIMUM`.
- Add direct boundary accuracy, certificate replay, cell/path traversal, and
  continuity tests.

### Phase 2 - Recipe and V2 evaluator

- Add recipe PODs, field-id/error transport, table-driven canonicalization, and
  deterministic curve/harmony resolution.
- Add the basis-aware 36-byte control snapshot and explicit `LEGACY_SINE`
  evaluator without changing legacy output.
- Implement unwrapped hue direction, local-gamut chroma, path-minimum chroma,
  explicit torsion, easing, loop layout, and direct/baked endpoint aliasing.
- Keep the legacy constructor untouched.

### Phase 3 - Morphing

- Make local-gamut `lerp()` interpolate relative chroma.
- Add global hue-phase alignment, achromatic hue inheritance, exact endpoint
  copies, and the restricted `PATH_MINIMUM` compatibility rule.
- Implement the allocation-free `PaletteCrossfade` dispatcher, rebake/view
  APIs, cancellation, and atomic target-policy handoff.
- Reconfirm both animation layouts and timeline capacity.

### Phase 4 - WASM and palette tool

- Add the versioned compile, inspect, generate, and transition-frame bridges.
- Update Daydream `tools/palettes.html` and the three companion palette modules
  according to section 14, while preserving Procedural and Legacy Profiles.
- Add reducer/adapter/canvas tests, the browser smoke test, and JS/WASM parity
  tests over the full option matrix.

### Phase 5 - Effect migration

- Migrate one low-risk baked consumer first.
- Review stills and transitions on host and device.
- Migrate remaining effects individually; do not batch-update visual goldens
  without inspection.

### Phase 6 - Cleanup

- Remove JS mirrors made obsolete by the V2 recipe bridge.
- Remove legacy APIs only after repository-wide call-site elimination.

## 18. Test requirements

### 18.1 Unit invariants

- Every curve resolves its exact three expected nodes.
- Base hue edits do not change non-hue resolved controls.
- Identical resolved recipes produce identical keys and LUT bytes; identical
  generation requests produce identical resolved recipes.
- Changing `variation_seed`, `sequence_index`, or a draw-site id independently
  changes the corresponding hashed draws; each is present in the hash input.
- Complementary `SHORTEST` resolves to segment deltas `+0.5,-0.5`; forced
  directions resolve half-turn ties with their specified sign and leave a
  zero wrapped delta at zero.
- Signed sweeps obey the exact `SHORTEST`/clockwise/counterclockwise interaction
  table; full-turn and multi-turn sweeps retain their signed winding.
- Harmony and custom loops resolve an explicit closing delta, and every segment
  reports the winding class given by the section 8.1 formula.
- Direct `get(1)` for single- and multi-turn loops is bit-identical to `get(0)`;
  a baked loop copies entry 0's RGB and alpha bits into entry 255.
- `MIRROR` is byte-symmetric around `t=0.5`.
- `LOOP` has identical first/last RGB. With cosine or smoothstep easing, its
  one-sided derivatives also match within LUT quantization tolerance.
- Constant-lightness recipes keep `L` constant before quantization.
- `LOCAL_GAMUT` tracks requested `q` within the boundary solver tolerance.
- `LOCAL_GAMUT` and `PATH_MINIMUM` both produce 62% of their respective
  boundary reference for `q=0.62, headroom=0.94`, and both clamp `q=1` to 94%.
- A `LOCAL_GAMUT` snapshot round-trips nonzero `q` and path hue exactly at
  both `L=0` and `L=1`; realized output is black/white without mutating them.
- `PATH_MINIMUM` keeps absolute `C` constant for a constant chroma curve.
- The committed interval certificate replays successfully and covers every
  closed LUT cell and radial interval `[0,floor]`; changing a coefficient,
  gamut constant, grid dimension, floor, or manifest hash invalidates it.
- Every point on a `PATH_MINIMUM` path is covered by a traversed certified LUT
  cell; stops and analytic lightness extrema are included even when they do not
  land on `i/255`.
- Torsioned recipes store untorsioned path hue, apply torsion exactly once,
  and derive or query all gamut quantities at the resulting final hue.
- `OKLAB_CARTESIAN` complementary paths approach a lower-chroma midpoint than
  their `OKLCH_ARC` equivalents.
- Vignette/falloff black stops remain exactly black.
- `falloff_start` accepts finite values strictly inside `(2/3,1)`, rejects both
  endpoints and non-finite values, and unequal starts are morph-incompatible.
- An OKLCH segment with exactly one `q=0` endpoint uses the chromatic endpoint's
  hue throughout; two `q=0` endpoints use zero. Injected black stops follow the
  same rule.
- `PATH_MINIMUM` rejects Cartesian, vignette, and falloff combinations.
- Every numeric field is tested at its finite boundaries, immediately outside
  them, and with NaN and both infinities. Native and WASM produce the same
  canonical recipe, adjustment masks, error code, and field id; failed calls
  leave output palettes and LUT buffers unchanged.
- Unknown enum values, excessive sweep/custom winding, invalid loop sweeps,
  negative ranges/chroma, and incompatible option combinations follow the
  exact table in section 7.1.
- Every pinned legacy constructor and `from_hsv_keys()` LUT remains byte-equal
  through the explicit `LEGACY_SINE` evaluator.

### 18.2 Exhaustive gamut census

Sweep:

- all 256 base hues;
- every harmony;
- representative constant/ascending/bell lightness curves;
- muted, balanced, vivid, and maximum chroma;
- every V2 domain and both color paths.

For non-black `LOCAL_GAMUT` samples using `OKLCH_ARC` and `headroom < 1`:

- direct linear RGB must be in gamut without invoking fallback mapping;
- no channel may contain a boundary plateau caused by chroma projection;
- any 8-bit zero caused solely by very low lightness is reported separately
  from gamut contact.

Maximum mode is expected to touch bounds and is excluded from the no-contact
gate. `OKLAB_CARTESIAN` is measured separately because an OKLab line between
two in-gamut keys is not guaranteed to remain in the sRGB gamut; fallback
mapping is legal there and its incidence is reported.

### 18.3 Morph tests

- Compatible wipes land exactly on target keys.
- Equal-`q` endpoints retain equal `q` through the wipe within tolerance.
- Wipes from `L=0` and `L=1` local-gamut keys retain their latent `q` while
  lightness moves into the chromatic range.
- A wipe from `q=0` blooms only in the target hue; a wipe to `q=0` retains the
  source hue until the exact target snapshot is copied. Both endpoints compare
  byte-equal to their respective input snapshots.
- Compatible custom palettes `[0,.25,.5]` and `[10,10.25,10.5]` produce no
  global spin; keys B and C preserve their coherently interpolated offsets from
  the selected chromatic anchor.
- Hue seams do not reverse direction or pop at half-turn crossings.
- Full-turn and multi-turn wipes preserve their winding count.
- Falloff palettes with unequal `falloff_start` values reject key morphing and
  reach the target exactly through baked output crossfade.
- `PATH_MINIMUM` palettes with bit-identical aligned geometry and
  `C_path_max` remain in gamut throughout a key wipe; any geometry difference
  rejects key morphing and dispatches a crossfade.
- Incompatible policies are rejected for key morphing before scheduling and
  dispatch directly to `PaletteCrossfade`.
- Crossfade endpoints are exact; completion copies the complete target policy
  into the stable subject address, cancellation restores the unchanged source,
  and pausing changes neither weight nor output.
- Crossfade stepping performs no arena allocation, does not call
  `bake_blend()`, and never reads either endpoint after its documented lifetime.
- Pause, first-step snapshot, and slow-fade resolution contracts remain intact.

### 18.4 WASM parity

- V2 bridge LUT bytes equal native engine samples at pinned recipes.
- All enum values and recipe schema versions are checked at the boundary.
- Diagnostics equal engine intermediate values within float tolerance.
- `generateV2` returns the same canonical recipe and status as the native
  generator for pinned requests and all 256 sequence indices.
- `bakeTransitionFrameV2` returns the same transition kind, incompatibility
  reason, and LUT bytes as native animation evaluation at pinned progress
  values for compatible and incompatible pairs.
- Mirror and loop endpoint handling is exact in the exported byte LUT.
- Loop entry 255 is a byte copy of entry 0 in both native and WASM output.
- Legacy `bakeLut` retains enum values `0..3`, and V2 domain values are accepted
  only by the versioned V2 bridge.

### 18.5 Performance and memory gates

Measure on Teensy, not only native host:

- A 256-entry V2 rebake must be no more than 1.25x the saturated legacy rebake
  cycle count unless the affected effect still clears its frame budget with
  documented margin.
- Sparse direct `get()` must be no more than 1.25x the saturated legacy path.
- Migrated pixel-hot effects must read a `BakedPalette`, not V2 directly.
- No migrated effect may regress its measured frame time by more than 2% from
  palette work.
- `ColorWipe` and `PaletteCrossfade` each stay inside the 112-byte device
  timeline slot.
- A 256-entry `rebake_crossfade` is measured separately and must be no more than
  2.25x a single V2 rebake; any caller must still clear its frame budget.
- V2 adds no persistent arena allocation, and a crossfade reuses the caller's
  existing display LUT without changing arena high-water usage.
- Record flash, ITCM, RAM1, and RAM2 deltas; any ITCM increase above 4 KiB is a
  stop-and-review condition, not an automatic acceptance.

### 18.6 Variety census

For every reference family, generate one full 256-step hue sequence and
report:

- unique baked-LUT hashes;
- lightness range, chroma range, and gamut-fraction range;
- fraction of samples classified achromatic;
- fallback gamut-map incidence;
- pairwise delta between consecutive palettes at the three resolved color keys.

The automatic base-hue sequence must visit all 256 8-bit hue values exactly
once before repeating. A family not explicitly intended to produce solid or
achromatic output must produce at least 240 unique LUT hashes over that cycle
and no solid-color LUT. These are degeneracy gates, not claims that a numeric
distance alone measures artistic quality. The same census produces a contact
sheet for visual review.

## 19. Acceptance criteria

V2 is ready for production migration when:

1. The default local-gamut mode does not use fallback gamut reduction in the
   exhaustive non-black census.
2. Constant-lightness and constant-chroma reference recipes meet their stated
   invariants.
3. The spectral loop has bit-identical direct and baked endpoints in native and
   WASM builds, including multi-turn sweeps.
4. `ColorWipe` and `PaletteCrossfade` remain layout-safe, deterministic, and
   visually continuous; crossfade completion and cancellation satisfy the
   subject-policy handoff contract without arena growth.
5. Native and WASM validation return identical canonical recipes, all three
   adjustment masks, errors, and field ids for the boundary matrix in section
   18.
6. The browser exposes independent hue/lightness/chroma/path controls and its
   diagnostic curves come from engine data.
7. Legacy constructors retain their pinned output until explicitly removed.
8. Device cycle, frame-time, arena, and ELF gates pass.
9. At least five reference families in section 13 are visually reviewed on the
   physical display, including one isolight sweep, one tonal monochrome, one
   neutral bridge, and one loop.

## 20. Recommended first implementation slice

The smallest slice that proves the design is:

1. Extract `gamut_max_chroma` without changing existing output.
2. Add `PaletteRecipe` with harmony/sweep hue, constant/ascending/bell
   lightness, constant/bell local-gamut chroma, straight/mirror/loop shapes,
   and OKLCH arc interpolation.
3. Preserve the 36-byte snapshot and existing `ColorWipe` size.
4. Add the V2 WASM bake plus `L`, `q`, and `C_max` diagnostics.
5. Demonstrate three recipes in the tool: balanced analogous, isolight
   spectral loop, and tonal monochrome.

Neutral bridges, path-minimum chroma, custom keys, and effect migrations follow
only after that slice passes gamut, performance, and symmetry gates.
