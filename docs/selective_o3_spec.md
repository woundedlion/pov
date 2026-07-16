# Selective -O3 regions for the Phantasm device build — spec

**Status: implemented 2026-07-15 — minimal region set** (§6 gate extension,
§4.1 macros, R1 sink, the fused scan drivers from R2, and a
DisplacementField-local draw-path promotion), landed `47090612..7284bb09`. The full
R1–R6 set is **measured-unaffordable** in the current tree — see the filled
ledger and the implementation-outcome notes in §7. Measured results:
`docs/profiles/shipping/`.

## 1. Goal

The shipping Phantasm image is `-Os` because the full 26-effect roster overflows
FlexRAM at `-O3` (`platformio.ini` §4.1 size relief valve). The 2026-07-14
on-device profile sweep (`docs/profiles/README.md`) measured what that costs:
renders run **1.14×–2.04× faster at -O3**, and the win concentrates almost
entirely in a handful of per-pixel inner loops (the `filter_blend` leaf alone
runs ~1.7–2.4× faster).

This spec adds a mechanism to compile **only those inner loops at `-O3`** while
the rest of the image stays `-Os`, recovering most of the frame-cadence wins for
a small, gated ITCM cost.

### Non-goals

- **No global `-O2` (or any global level change).** Rejected by the owner.
- **No LTO**, no linker-script changes, no per-TU build splits (the device build
  is effectively a single TU — the converted `.ino.cpp` plus two small
  `core/engine/*.cpp` files that include no render headers).
- **No new FlexRAM bank for ITCM.** See §3 — this is a hard wall, not a budget.
- **No flash-resident hot code.** Moving `-O3` overflow to XIP flash would trade
  deterministic POV column timing for cache-miss jitter; out of scope.
- **No change to any non-Phantasm image**: `holosphere` (already `-O3`), host
  tests (clang), and WASM (emscripten) must be bit-for-bit unaffected.

## 2. Measured motivation (what to optimize for)

Frame cadence quantizes to 62.5 ms display windows, so a speedup only matters
when it crosses a window boundary. From `docs/profiles/README.md`, the effects
that gain (or partially gain) a cadence tier at -O3, and the hot scope that must
speed up to buy it:

| Effect | Cadence -Os → -O3 | Hot scope (profile counter) |
|---|---|---|
| MeshFeedback | 8 → 16 fps | feedback flush (`mf_feedback_flush`) |
| Flyby | 8↔16 → steady 16 fps | shader (`fly_shader_draw`) |
| Voronoi | 8 → 16 fps | per-pixel KD walk (`vo_shade`) |
| DreamBalls | 5.3 → 8 fps (P0) | wireframe raster (`db_mesh_plot`) |
| FlowField | 5 → 8 fps (saturated) | particle raster |
| HopfFibration | 8↔16 jitter → steady 16 | trail raster |
| IslamicStars | 8 → 16 fps (mid shapes) | per-face SDF (`is_mesh_scan`) |
| HankinSolids | crossfade tier gain | per-face SDF (`hk_mesh_scan`) |
| RingSpin | 8 → **16 fps locked** (fused tip, 2026-07-15) | fused ring-group raster (`rs_ring_scan`) |
| DisplacementField | 8 → 16 fps, whole cycle (fused tip, 2026-07-15) | fused stack scan (`df_fused_scan`) |

Effects already at 16 fps at -Os (Comets, PetalFlow, RingShower, Thrusters,
MobiusGrid, ChaoticStrings, GnomonicStars, SphericalHarmonics, Liquid2D) gain
**nothing user-visible** from -O3 — their savings fall into `*_buffer_wait`
idle. Raymarch, MindSplatter, Dynamo, BZ/GS stay 8 fps even at full -O3. **Do
not spend ITCM on code that only serves these effects.** (They may still benefit
incidentally when they share a region with a tier-crossing effect; that is fine —
it just isn't a reason to add a region.)

The optimization objective is therefore: **maximize cadence-tier crossings per
ITCM byte**, in the table's order of value.

## 3. Memory model and the hard budget

Teensy 4.0 FlexRAM = 512 KiB in sixteen 32 KiB banks, partitioned at boot into
ITCM (code, rounded **up** to whole banks) + DTCM (variables + stack). The
calibrated shipping figures (`tools/teensy_budgets.json`, phantasm):

```
RAM1 used 509,312 B, free-for-locals 14,976 B, stack floor 12,288 B
⇒ ITCM code 146,456 B → 5 banks = 163,840 B → intra-bank padding 17,384 B
⇒ DTCM 11 banks = 360,448 B: variables 345,472 B + free 14,976 B
```

Two consequences the implementer must treat as invariants:

1. **A 6th ITCM bank is a hard failure, not a squeeze.** It would shrink DTCM to
   10 banks = 327,680 B < the 345,472 B of variables — the image cannot even
   link/boot, independent of the stack floor. The **only** growth room is the
   intra-bank padding (~17.4 KB at calibration).
2. **The existing gate cannot see this growth.** `teensy_size` reports
   `RAM1: variables + code + padding`; code growing into padding moves bytes
   between two components of the *same sum*, so RAM1 `used` and `free` are both
   unchanged until the bank cliff — where the build breaks all at once. §6 adds
   a component-level ceiling to make the growth visible and ratcheted.

One property makes growth-into-padding strictly safe rather than a gamble:
ITCM is tightly-coupled memory — single-cycle, not cached — so -O3 code added
to a bank cannot evict, pressure, or slow down any other code. This is why the
"incidental ride-along" effects in §2/§5 are free, and why the only cost
dimension of a region is bytes; there is no I-cache-pressure counterargument.

Numbers above are calibration-time; the working tree moves. **Phase 0 (§7)
re-measures them on the actual base commit before any region lands.**

Growth reserve rule: after the final region set, ITCM padding must remain
**≥ 4,096 B** (headroom for ordinary future code growth before anyone has to
think about banks again). This bounds total -O3 growth to roughly 13 KB at
calibration figures.

FLASH is not a constraint (~1.8 MB free); every ITCM function also has a flash
load image, so expect flash to grow by about the same amount plus any literal
pools — ignore it except for the existing gate ceiling.

## 4. Mechanism: `HS_O3` region macros

### 4.1 Definition

Add to `core/engine/platform.h`, directly below the `HS_COLD` block
(currently ends near line 1019), with a comment in the same style:

```c
// ---------------------------------------------------------------------------
// HS_O3_BEGIN / HS_O3_END: compile the enclosed function definitions at -O3 on
// the -Os device image (selective hot-loop optimization; docs/selective_o3_spec.md).
// Active only for device GCC building at -Os (__OPTIMIZE_SIZE__): the holosphere
// -O3 image, host clang, and WASM see no-ops, so those builds are byte-identical.
// The fast-math flags are restated because GCC 11's optimize pragma rebuilds
// optimization flags from defaults, dropping command-line -ffast-math /
// -fno-finite-math-only for the region (fixed in GCC 12; harmless to restate).
// HS_O3_FN is the single-function fallback for definitions a region cannot wrap.
// ---------------------------------------------------------------------------
#if defined(ARDUINO) && defined(__GNUC__) && !defined(__clang__) && \
    defined(__OPTIMIZE_SIZE__)
#define HS_O3_BEGIN                                                     \
  _Pragma("GCC push_options")                                           \
  _Pragma("GCC optimize(\"O3\", \"fast-math\", \"no-finite-math-only\")")
#define HS_O3_END _Pragma("GCC pop_options")
#define HS_O3_FN __attribute__((optimize("O3", "fast-math", "no-finite-math-only")))
#else
#define HS_O3_BEGIN
#define HS_O3_END
#define HS_O3_FN
#endif
```

Activation matrix (must hold; §8 verifies):

| Build | Compiler | Level | Macros |
|---|---|---|---|
| `phantasm`, `phantasm8`, `profile` | arm-gcc 11.3.1 | `-Os` | **active** |
| `holosphere`, `holosphere_dma`, `profile_o3` | arm-gcc 11.3.1 | `-O3` | no-op (`__OPTIMIZE_SIZE__` undefined) |
| native tests | clang | any | no-op |
| WASM release/debug | emscripten clang | any | no-op |

(`defined(ARDUINO)` makes the matrix true by construction — a hypothetical
host-GCC `-Os` build can never activate regions — mirroring how the `HS_COLD`
block relies on the platform fork.)

Two useful properties fall out: `just profile <Effect>` (the `-Os` `profile`
env) measures the **shipping** selective-O3 configuration, and `profile_o3`
remains an untouched **global-O3 reference** to compare against.

Form choice: prefer the **region** form for clusters of adjacent definitions
that must share one option node (matching options inline into each other
freely) and for anything containing lambdas or member templates — exactly
where per-function attribute application is historically flaky. Prefer
`HS_O3_FN` for a single ordinary function: it cannot be left unmatched and
cannot accidentally span an `#include`. Both produce the same option node; the
choice is ergonomics and failure surface, not semantics.

### 4.2 GCC 11.3.1 semantics the implementer must respect

These are load-bearing; getting any one wrong silently produces a slower or
semantically different image.

1. **Flag reset (why fast-math is restated).** Before GCC 12, `#pragma GCC
   optimize` / `__attribute__((optimize))` rebuild the optimization-flag set
   from the *defaults* of the requested level, not from the command line. A bare
   `optimize("O3")` region would drop `-ffast-math -fno-finite-math-only` inside
   the region. Worse, dropping only `-fno-finite-math-only` while keeping
   `fast-math` would fold `std::isfinite()` to constant-true in the render
   sink's NaN guard (`core/render/filter.h`, see the platformio.ini §4.1 flag
   comment). The macro therefore restates both, exactly as written above. Do
   not add or remove flags without re-running the §8 checks — **and do not
   reorder them**: the list is processed like a command line, and
   `no-finite-math-only` must follow `fast-math` (which implies
   finite-math-only) or the isfinite guard folds anyway. One further named
   risk: `fast-math` is a driver-umbrella token, and GCC's optimize pragma may
   *ignore* rather than expand it — then the region silently runs without
   fast-math: safe (the guard stays live) but slower, reading as an "-Os wall".
   The §8 codegen checks are the detector for both directions — the
   differential / `profile_o3` comparison for a dropped umbrella, the isfinite
   probe for a folded guard. If the umbrella proves ignored, restate its
   constituent sub-flags instead (`no-math-errno`,
   `unsafe-math-optimizations`, `no-signed-zeros`, `no-trapping-math`).
2. **Inlining across an optimization boundary is restricted, not forbidden.**
   GCC's mismatch refusal keys on the callee *carrying its own* optimize
   options or semantically incompatible FP flags: a plain, attribute-free
   `inline` helper in the same TU is generally still inlinable into a region
   caller, and the inlined body is re-optimized under the caller's -O3 options.
   (This is a second reason the fast-math restatement is load-bearing: matching
   FP flags keep the inline compatibility check passing.) `always_inline`
   helpers always cross and adopt the caller's options. Consequence: regions
   can likely stay **narrow loop drivers**, with the shared plain-inline
   helpers (`blend_alpha` `core/color/color.h:411`, `quintic_kernel`
   `core/math/3dmath.h:72`, `fast_wrap` `core/engine/util.h:86`, the Canvas
   accessors) folding in for free — but this must be confirmed by the Phase 1.5
   probe (§7) before sizing any region. If the probe fails, each region must
   transitively cover its helpers (contingency `HS_O3_FN` in color.h /
   3dmath.h / util.h / sdf.h) and per-region ITCM estimates roughly double.
   A hot function *called from* -Os code is fine either way (that boundary
   call is per-draw, not per-pixel).
3. **The reverse inlining direction silently evaporates a region.** An inline
   region function that is itself inlined into an -Os *caller* adopts the
   caller's options; the pragma protects only the out-of-line copy. The big
   rasterizer bodies are far over -Os inline limits, so this is unlikely in
   practice — but the §8 "pragma applied" checks (map diff, disassembly) both
   inspect the out-of-line symbol and would pass even if a hot call site had
   inlined an -Os copy. The sentinel profile delta is the real detector for
   this failure mode; do not trust the disassembly check alone. The
   complementary outcome at the same boundary is the **seam**: once a function
   carries region options, -Os callers that previously inlined it may now call
   it out-of-line (per #2's mismatch rule). The out-of-line body is -O3, but
   the call overhead and the lost cross-function optimization mean a region
   measured in isolation can undersell the uniform-O3 ceiling — seams heal
   when the adjacent caller region lands with *matching* options. §7's
   decision rule accounts for this (R1 is not reverted on its solo result).
4. **Templates take the options at their definition site.** Wrapping a template
   definition in a region covers **all instantiations**, including per-effect
   typed-pipeline instantiations made elsewhere. This is how one region serves
   many effects — and also how one region multiplies its size cost (§5 flags
   which regions instantiate per-effect).
5. **Region hygiene.** Between `HS_O3_BEGIN` and `HS_O3_END`: complete function
   /struct definitions only. No `#include`. No functions marked `HS_COLD` /
   `HS_COLD_MEMBER` (contradictory placement/optimization intent). Every
   `HS_O3_BEGIN` must have a matching `HS_O3_END` in the same file; an unmatched
   push corrupts option state for the rest of the TU (which, here, is the whole
   image).
6. **Verify application empirically.** GCC has historically had holes in pragma
   application (lambdas, templates in some versions). The per-region procedure
   (§7) treats "ITCM code delta ≈ 0 **and** no profile delta" as "the pragma did
   not take": first check the map for the expected symbols growing, then fall
   back to `HS_O3_FN` on the specific function(s) instead of the region form.
7. **No ODR hazard in practice** — the device image is one TU, and every other
   build sees no-op macros — but keep regions in headers self-contained anyway
   (property 4) so two device TUs could never see different option states for
   the same inline function. (`core/engine/memory.cpp` / `reaction_graph.cpp`
   include no render headers today.)

## 5. Region inventory

Anchors are from the 2026-07-14 tree; re-locate by symbol name if lines have
drifted. Apply **in this order** (shared-instantiation, high-fan-out regions
first), measuring each per §7 before starting the next.

### R1 — blend leaf + AntiAlias splat (`core/render/filter.h`)

The single highest-value region: the profiles show this leaf ~halving at -O3
(1.7×–2.4×). The terminal sink is genuinely **one instantiation per `<W,H>`**
serving every effect that plots through a pipeline (its out-of-line copy is
pinned by `PipelineRef` taking its address — which is also what makes the
Phase 1.5 objdump possible). The AntiAlias half is *not* one instantiation —
see its bullet.

- `Pipeline<W,H>` terminal specialization (`filter.h:100`), specifically its
  `plot(Canvas&, int x, int y, const Pixel&, float, float alpha)` — the
  `HS_PROFILE(filter_blend)` scope (`:132–140`) and the `blend_alpha` write.
  `blend_alpha` is a plain-`inline` functor factory at `core/color/color.h:411`
  (not `always_inline`); per §4.2 #2 it is expected to inline into the O3 sink
  and re-optimize there — this is exactly what the Phase 1.5 probe confirms.
  Contingency if it does not: `HS_O3_FN` on `blend_alpha` in color.h, not a
  wider R1.
- `Screen::AntiAlias<W,H>::plot` (`filter.h:965/981–1028`) — the 4-tap quintic
  splat feeding up to four `filter_blend` writes per sample. Two caveats,
  both measured-not-assumed: it is a member template on the forwarding
  `PassFnT` callback, so it instantiates **once per distinct downstream chain
  (per-effect)**, multiplying like R3; and its doc comment (`filter.h:978`)
  says the forwarding design exists precisely so the taps inline into the
  caller — wrapping it makes it option-carrying, which may stop that inlining
  (a §4.2 #3 seam) until the caller loop joins a matching-options region (R3).

Serves: every tier-crossing effect except Flyby and Voronoi (which bypass the
pipeline sink). **Do not** wrap the recursive `Pipeline<W,H,Head,Tail...>` node
(`filter.h:281–309`) in R1 — it instantiates per-effect; its per-pixel body is
a thin router that the profiles do not flag. Revisit only with measurements.

Because of the §4.2 #3 seams, R1's solo sentinel result is a **lower bound**:
do not revert R1 on a weak solo number — record it in the ledger and defer its
final keep/revert verdict until after R3 lands and the seams heal.

### R2 — SDF scan kernel (`core/render/scan.h`)

Shared per-`<W,H,ComputeUVs>` by design (type-erased shader, `PipelineRef`
default — see the comment at `scan.h:46–48`); typed-pipeline callers add
instantiations (RingSpin, DistortedRing, ShapeShifter — watch the map).

- `Scan::process_pixel` (`scan.h:50–118`) — per-pixel SDF eval + AA coverage +
  shader + plot.
- `scan_region` (`scan.h:143`) and `Scan::rasterize` (`scan.h:356–388`) — the
  loop drivers, wrapped so the pixel loop and `process_pixel` inline together.
- **The fused scan drivers (landed after this spec's first draft; both are
  self-contained per-pixel loops that bypass `process_pixel`):**
  `Scan::DistortedRingStack::draw` (DisplacementField's same-axis stack scan;
  its per-pixel frame + candidate loop calls
  `SDF::DistortedRing::distance_from_frame` / `polyline_distance`) and
  `Scan::RingGroup::draw` (RingSpin's per-trail-frame group scan; its pixel
  lambda calls `SDF::Ring::distance<false>` per slot — an out-of-line call at
  -Os, which is exactly the per-blend chain the 2026-07-15 profile measured
  at ~960 cyc/blend). Their pixel lambdas are defined inside the wrapped
  functions, so they inherit the region options (§4.2 #4); the distance
  callees in sdf.h must be confirmed inlined per the R2 post-land check
  below. These two drivers are the measured highest-value additions: the
  2026-07-15 fused-tip profiles show they carry RingSpin's full 8 → 16
  crossing and DisplacementField's whole-cycle crossing at -O3.

Serves: IslamicStars + HankinSolids (`Scan::Mesh::draw` at `scan.h:787–852`
calls `Scan::rasterize`; its per-face setup loop may join the region if the map
shows it matters), RingSpin (`RingGroup`), DistortedRing, DisplacementField
(`DistortedRingStack`), ShapeShifter.

The per-pixel arithmetic R2 exists to speed up lives in `core/render/sdf.h`,
reached as `shape.distance<ComputeUVs>()` from `process_pixel` — a `distance()`
implementation left at -Os is precisely the "-Os wall inside a loop" §9 warns
about and would gut the region. Expected coverage, per shape family:

- **`SDF::Face`** (IslamicStars/HankinSolids): its internals are saturated with
  `always_inline` (`sdf.h:2149–2578`), which crosses the optimization boundary
  unconditionally and adopts the region's -O3 — covered by construction.
- **Ring / knot-polyline SDFs** (RingSpin, DistortedRing, DisplacementField):
  plain-inline `distance()` members (`sdf.h:627/830`; `polyline_distance`
  `:898`) — expected to inline per §4.2 #2 / Phase 1.5.
- **Composite shapes** (`sdf.h:1144/1282/1492`): same expectation.

After R2 lands, confirm via map/objdump that no `distance()` symbol is called
out-of-line from the region; any that is gets its own region (or `HS_O3_FN`)
in sdf.h around that shape's `distance()`.

### R3 — plot rasterizer (`core/render/plot.h`)

- `Plot::rasterize` (`plot.h:680–928`) including its `process_segment` lambda
  (`:746`) — lambdas defined inside a region function inherit its options.
- `rasterize_planar_strategy` (`plot.h:182`) and `rasterize_geodesic_strategy`
  (`plot.h:294`).
- `Plot::Mesh::draw_edge` (`plot.h:2241–2267`) — thin, but it must not become
  an -Os wall between mesh draws and `rasterize`.

Instantiation risk: `PipelineT` defaults to shared `PipelineRef`, but effects
passing typed `filters` instantiate per-effect — this region has the largest
multiplication potential. If the map shows unacceptable growth, the fallback is
`HS_O3_FN` on the two strategy leaves only (they are the inner loops; the outer
`rasterize` body is per-segment, not per-pixel).

Serves: DreamBalls, HopfFibration, FlowField, MeshFeedback's wireframe pass,
plus most 16-fps effects incidentally.

### R4 — shader compositors (`core/render/scan.h`)

- `Scan::Shader::draw`, closure-shader variant (`scan.h:954`) — instantiates
  per-effect (`ShaderFn` is the effect's lambda: Flyby `Flyby.h:110–125`,
  Liquid2D `Liquid2D.h:132`). Wrapping the definition covers both; the effect
  lambdas inline into the O3 loop and inherit it.
- `Scan::Shader::draw`, split fragment/vertex variant (`scan.h:1017`) — shared
  per `<W,H,SAMPLES>` (type-erased shaders; BZ/GS). BZ/GS cross no tier.
  **Default: leave this overload at -Os**; promote it only if the map diff
  shows it sharing sections with the closure variant or costing ~nothing.

Serves: Flyby (tier crosser). Liquid2D rides along (1.14×, no tier — fine,
it's the same instantiation site).

### R5 — feedback/trails flush (`core/render/filter.h`)

MeshFeedback is the single biggest tier win (2.04×, 8 → 16 fps); Dynamo shares
the machinery (no tier, rides along).

- `Pipeline<...>::flush` chain (`filter.h:405–426`), `Screen::Trails::flush`
  (`filter.h:1100`), and the terminal `Pixel::Feedback::flush` (`filter.h:1283`)
  whose per-pixel `composite` loop (`filter.h:1403–1490`) dispatches the blend —
  the actual hot body the profiles count. The effect-side color lambdas
  (`MeshFeedback.h:154–160`, `Dynamo.h:277–282`) inline into the O3 flush loop.
- The flush machinery is per-effect-instantiated (typed filter lists) but only
  a few effects carry flush-bearing filters; the map diff will show exactly
  which instantiations grow.

### R6 — Voronoi KD walk (effect-local)

Voronoi bypasses the pipeline entirely (writes `canvas(x,y)` directly,
`Voronoi.h:244`), so R1–R5 do not touch it.

- The `vo_shade` double loop (`effects/Voronoi.h:199–246`).
- `KDTree::nearest` (`core/mesh/spatial.h:108`) — must join the region (or get
  `HS_O3_FN`) or the per-pixel `tree.nearest(p,2)` call at `Voronoi.h:235`
  becomes a non-inlined -Os call and the region buys little.

### Explicitly out (do not add without new profile evidence)

BZ/GS `step_physics` sim kernels, Raymarch's ray-march (`Scan::Volume::draw`,
1.23×, no tier), MindSplatter, GnomonicStars' star transform, and all
effect-local code of the always-16-fps effects. If budget remains after R1–R6,
prefer **deepening** existing regions (e.g. `Scan::Mesh::draw`'s face-setup
loop) over new effect-local ones.

## 6. Size-gate extension (prerequisite, land first)

Make ITCM code growth visible and ratcheted before any region lands.

1. **Schema** (`tools/teensy_budgets.json`): allow an optional per-component
   ceiling inside a region entry:

   ```jsonc
   "ram1": {
     "max_bytes": 512000,
     "free_min_bytes": 12288,
     "components": { "code": { "max_bytes": <calibrated> } }
   }
   ```

2. **Gate** (`tools/teensy_gate.py`): `parse_teensy_size` already captures the
   per-region `components` dict; extend `evaluate()` to check
   `regions.<r>.components.<name>.max_bytes` against
   `sizes[r]["components"][name]`, with a new violation code
   `component-over-budget`. A configured component **missing** from the parsed
   output is a hard `component-missing` violation (same fail-loud rule as
   regions/symbols — a renamed teensy_size field must not silently disable the
   ceiling).
3. **Tests** (`tools/teensy_gate_tests/`, auto-discovered by
   `just teensy-gate-test` and CI): fixture + cases for pass, over-ceiling, and
   missing-component. Follow the existing golden/deliberately-broken fixture
   pattern.
4. **Calibration and the per-commit ratchet**: set the phantasm
   `ram1.components.code.max_bytes` to `Phase-0 measured code + 2,048` when §6
   lands (a pure ratchet proving the gate bites before any region exists).
   **Every kept region's commit then bumps the ceiling to its own measured
   code + 2,048** — the budgets file's standard "raising a ceiling is a
   reviewed one-line edit landed with the change that needs it" rule. Without
   the per-commit bump, the first region commit turns the gate red. The final
   value must also satisfy the §3 reserve rule, expressed bank-relative so it
   survives baseline drift: `bank_ceil(code) − code ≥ 4,096`, with the ITCM
   bank count unchanged from Phase 0. Holosphere's budget entry is unchanged
   (no components key — the feature is opt-in per target).
5. Update `docs/teensy_ci_gate_spec.md` (§8 budgets schema) to document the
   `components` key.

## 7. Implementation procedure

The HS_PROFILE instrumentation and the `profile`/`profile_o3` envs are now
landed (`9c994a1e`), and the scopes compile to nothing without
`HS_PROFILE_ENABLE` — the shipping `phantasm` image is unaffected, so size
measurements and profiling runs both work from clean commits. Land via the
usual worktree + FF integrator workflow, one region per commit.

**Phase 0 — baseline.**
1. `pio run -e phantasm` on the base commit; record the `teensy_size` RAM1 line
   (variables / code / padding / free) and keep `firmware.map`
   (`.pio/build/phantasm/firmware.map`, emitted by `tools/teensy_map.py`).
2. Recompute the §3 bank math from the measured line. If padding < 6,144 B —
   the 4,096 B reserve plus ~2 KB, the smallest region set worth landing —
   stop and renegotiate the reserve rule with the owner before proceeding.
3. Land §6 (gate extension + first calibration).
4. Capture -Os sentinel profiles on the base commit if the
   sentinel-baselines figures no longer match the tree (`just profile
   MeshFeedback|Flyby|Voronoi|DreamBalls|IslamicStars|RingSpin` — RingSpin
   included: it is R1's primary sentinel).

**Phase 1 — macros.** Add the §4.1 block to `platform.h`. Verify the
activation matrix cheaply: `pio run -e phantasm` and `-e holosphere` produce
byte-identical images vs Phase 0 (macros defined but unused), and `just test` +
`just build` still pass (no-op expansion on clang).

**Phase 1.5 — inlining probe (determines the sizing strategy for all six
regions).** Wrap a single minimal region — the R1 sink alone
(`Pipeline<W,H>::plot`) — build phantasm, and objdump the out-of-line sink
symbol (locate it via the map): confirm `blend_alpha` and the clip/saturation
helpers **inlined into it with -O3 codegen** (unrolling, FMA contraction), not
`bl` calls back into -Os code. If plain-inline helpers cross the boundary
(expected per §4.2 #2), keep every §5 region as drawn — narrow loop drivers.
If they do not, each region must transitively include its helpers (`HS_O3_FN`
on `blend_alpha`, `quintic_kernel`, `fast_wrap`, and the per-shape
`distance()` members in sdf.h), the per-region ITCM estimates roughly double,
and the plan must be re-checked against the §3 budget before Phase 2 starts.

**Phase 2 — regions, one at a time, in §5 order.** Per region:

1. Apply the region; `pio run -e phantasm`.
2. **Size**: record the RAM1 `code`/`padding` delta; diff the map against the
   baseline to confirm growth is the expected symbols (and to catch per-effect
   instantiation blowups). If code delta ≈ 0, suspect pragma non-application
   (§4.2 #6) before concluding the region is free.
3. **Speed**: `just profile <sentinel>` for the region's tier-crossing
   sentinel(s) (R1: RingSpin or DreamBalls; R2: IslamicStars; R3: DreamBalls;
   R4: Flyby; R5: MeshFeedback; R6: Voronoi). Compare the hot-scope figure and
   cadence against the sentinel-baselines table below (distilled from the
   `docs/profiles/shipping|O3/` reports).
4. **Decide**: keep if it buys a cadence tier, or closes **≥ half the remaining
   gap** between the -Os baseline and the -O3 ceiling figure for a
   tier-crossing sentinel's hot scope (see the baselines table below), at
   acceptable cost. The half-gap rule does not contradict §2's
   tiers-are-all-that-matters thesis — it exists because regions **stack**
   within one effect's frame (IslamicStars needs R1 + R2; neither alone may
   cross): a partial win on a *shared* region (R1–R3) is retained because a
   later region completes the crossing. It does **not** apply to effect-local
   regions (R6) — those either cross their tier or revert. Revert if the map
   shows multiplication out of proportion to the win (exception: R1's solo
   verdict is deferred past R3, §5). Before reverting a real-but-too-expensive
   win, try the -O2 variant (§9). Record every region in the ledger below —
   including reverted ones, with the measured reason — so nobody re-attempts a
   measured-dead region.
5. Stop when the §3 reserve rule would be violated, or when remaining regions
   serve no tier crossing.

Ordering tension, acknowledged: the three headline 8→16 effects in §10 depend
on R4/R5/R6 — the *last* regions applied. Shared-first ordering is still
correct (R1–R3 are the cheapest bytes per effect served and de-risk the budget
early), but if the budget exhausts before R4–R6, re-run the keep/revert
decision **globally**: a kept shared region that turned out to serve no
remaining tier crossing can be traded out for a headline region.

**Phase 3 — final calibration + full gates.** Re-set the component ceiling
(§6.4); `just teensy-size` (all four envs) green; final sentinel sweep.

### Ledger (filled 2026-07-15)

Phase 0 tripped the §7 stop: the DisplacementField-era landings had eaten the
calibration padding to **2,824 B** (code 161,016). Unblocked per owner
decision by (a) a cold-code ITCM eviction sweep (`2c2470b2`, −5,600 B) and
(b) the owner's Phantasm playlist trim of Dynamo + Thrusters (`2d21f8a4`,
−5,936 B), giving 14,360 B of padding at the pre-region tip (code 149,480).

| Region | ITCM code Δ (B) | Padding left (B) | Sentinel result | Kept? |
|---|---|---|---|---|
| R1 filter.h blend sink (sink only, no AntiAlias) | +720 | 13,640 | judged with R2 (seam rule) | ✅ `683d0eda` |
| R2 fused drivers only (`RingGroup`, `DistortedRingStack`) + sdf.h `HS_O3_FN` (`distance_from_frame`, `polyline_distance`) | +2,784 | 10,856 | RingSpin scan 34.6 ms vs 33.3 ceiling (+4%), **16 fps locked all windows**, blend 84 cyc | ✅ `0f7b0616` |
| DF-local: `draw_rings` `HS_O3_FN` (bake + hue prep + shader-lambda seam) | +3,488 | 7,368 | dwell scan 43–48 ms vs 44–46 ceiling; most of cycle locks 16 fps | ✅ `e98f4651` |
| DF-local: hue-table members `HS_O3_FN` | +736 | 6,632 | measured dead — recapture identical within noise (the cost is the still--Os OKLab chain in color.h, not the member bodies) | ❌ reverted `7284bb09` |
| R1 AntiAlias + R2 shared kernel (`process_pixel`/`scan_region`/`rasterize`) | +27,872 over the minimal set (probe, full R1+R2 = +31,360) | n/a — crosses the 6th-bank hard wall | not attempted on device | ❌ measured-unaffordable |
| R1 AntiAlias splat (region form; heals the R1 sink seam — the -Os splat could no longer inline the option-carrying sink) | +1,888 | 12,824 | judged with R3 (seam rule) | ✅ 2026-07-15 |
| vector_to_pixel `HS_O3_FN` (geometry.h; shared per `<W,H>`) | +336 | 12,488 | judged with R3 | ✅ 2026-07-15 |
| R3 geodesic strategy region (`rasterize_geodesic_strategy` + sampler lambda) | +3,296 | 9,192 | HopfFibration median per-call 220 → ~155 µs (at O3 ceiling); peak windows still spill | ✅ 2026-07-15 |
| R3 `Plot::rasterize` region (whole sim+draw loop; strategies + AntiAlias + sink fuse in) | **−6,656** (one shared -O3 copy per pipeline type replaces per-caller -Os inlining) | 15,848 | HopfFibration `hf_trail_raster` 220 → 144 µs/call median, 382 → ~210 µs peak-window — at/below the global-O3 ceiling (152/239); residual 2-window frames in the heaviest fold/twist phases are workload peaks the 68 s O3 reference pass never sampled, not an -Os wall. DreamBalls/MindSplatter/PetalFlow/MobiusGrid ride the same instantiations | ✅ 2026-07-15 |
| R4 shader compositors: closure `Shader::draw` `HS_O3_FN` + per-pixel callee chase (`stereo_noise_warp`, `SingleOpenSimplex2`, `hue_rotate_rgb`, `oklab_to_linear_rgb`, `linear_rgb_in_gamut`) | +3,104 | 12,680 (post mesh-scan driver + thunk flash-routing) | Flyby worst preset 50.9 ms vs 43.6 ceiling, **16 fps locked over the full cycle** (was 83.7 ms, 8 fps on 4 of 5 presets); loop+lambda alone and +noise-path were each measured ~flat — the OKLab hue chain was the cost | ✅ 2026-07-15 |
| MindSplatter wrapper: `Plot::ParticleSystem` region (both `draw` overloads — per-trail tween/cull/dispatch loop) | +3,312 | — | measured dead 2026-07-16 — per-preset scan identical within noise (worst 108.6 → 109.0 ms) | ❌ reverted, not landed |
| MindSplatter effect-local: `draw_particles` `HS_O3_FN` (mobius/hole/palette shader lambdas) | +1,184 | — | measured ~dead 2026-07-16 — uniform −1.2 % (worst 108.6 → 107.3 ms), no cadence change | ❌ reverted, not landed |
| `Plot::gate_trail_edges` region (hoisted per-trail clip gate: shared per-point rows/columns, whole-trail coarse reject, bits feed `rasterize`'s cull) + HopfFibration `gate_trails` wiring | +1,648 | 10,840 | HopfFibration: replaces the per-edge in-place gate; fully-invisible trails skip stage+rasterize whole | ✅ 2026-07-16 |
| R5 feedback flush | not measured | | | deferred — no budget |
| R6 Voronoi KD | not measured | | | deferred — no budget |

**MindSplatter's ship→ceiling ratio (1.15–1.49× per preset) is not
recoverable by selective -O3 — do not re-attempt.** With the R3 rasterize
chain, the `ParticleSystem` wrapper, and the effect shader lambdas all
measured, essentially no `-Os` code remains inside `msp_particle_scan`; the
residual gap is dominated by a cross-config workload artifact: both configs
build with `-ffast-math -fno-finite-math-only` (platformio.ini §4.1), but
`-O3` exercises the reassociation/FMA license far more aggressively than
`-Os`, so the physics produces different float trajectories — hence different
coverage — under identical flags. Not reachable compiler headroom; the
cadence lever is coverage (pool/trail size), not the compiler.

FLASH grew ~3–7 KB total across the kept commits — far under the gate ceiling
(not tracked per commit).

R5 and R6 are **deferred, not dead**: each needs ITCM headroom that the owner can
mint by further Phantasm playlist trims (precedent: `2d21f8a4` freed ~5.9 KB
by dropping two effects; `f6c2ff7c` earlier removed three whole features).
Budget the next region set from a fresh Phase-0 measurement, not this table's
padding column. Note the probe economics: the shared R2 kernel multiplied by
per-shape/per-pipeline instantiation to ~28 KB — narrow per-driver regions
plus targeted `HS_O3_FN` on the out-of-line hot callees bought the same
sentinel wins for a tenth of the bytes. §4.2 #2's inline expectation held for
most helpers (blend chain, quintic in the sink) but NOT universally: GCC 11
refused several plain-inline -Os callees inside regions
(`ClipRegion::render_y_end`, `polyline_distance`, the OKLab chain,
`quintic_kernel`/`wrap_t` at some sites) — expect to chase per-callee
`HS_O3_FN` from the disassembly, and judge by the sentinel counter, when
landing the next region.

Codegen checks (§8): (a) PASS — a throwaway bare-`optimize("O3")` build
produces different sink codegen, so the restated fast-math tokens reach the
region. (b) resolved-no-symbol — under global -O3 the sink has no out-of-line
copy (fully inlined), so equivalence was judged by the sentinel deltas.
(c) N/A — no live `isfinite` sits inside the landed regions (the Feedback
flush guard is R5, not landed); the restatement ordering is still verified by
(a). Holosphere: symbol inventory and all region totals byte-identical across
the region commits; the image hash differs only via shifted `__LINE__`
constants in always-on check macros (the region markers add source lines).

### Sentinel baselines (device captures; committed reports under `docs/profiles/{Os,O3}/`)

The full reports are now committed (one current report per effect per level);
the acceptance-relevant figures are distilled here as the durable referent.
Hot-scope render time per frame, shipping `-Os` vs the global `-O3` ceiling
(the full reports carry the window/cadence detail). RingSpin and
DisplacementField figures are from their 2026-07-15 **fused-tip** recaptures
(`56d8c854` / `7d50b672`); the rest are 2026-07-14:

| Sentinel | Hot scope | -Os | -O3 ceiling | Cadence -Os → -O3 |
|---|---|---|---|---|
| MeshFeedback | `mf_feedback_flush` (light morph) | 97.9 ms | 48.0 ms | 8 → 16 fps |
| Flyby | `fly_shader_draw` (expensive regime) | 77.0 ms | 42.4 ms | 8↔16 → 16 fps |
| Voronoi | `vo_shade` | 76.9 ms | 57.2 ms | 8 → 16 fps |
| DreamBalls | `db_mesh_plot` (P0, 18 copies) | 145.3 ms | 91.4 ms | 5.3 → 8 fps |
| IslamicStars | `scan_mesh_raster` | 60.3 ms | 43.4 ms | 8 → 16 fps (mid shapes) |
| RingSpin | `rs_ring_scan` (fused group, sustained) | 66.8 ms | 51.6 ms | ~9.4 fps → **16 fps locked** |
| DisplacementField | `df_fused_scan` (NOISE dwell) | 59.1 ms | 46.0 ms | 8 → 16 fps, whole cycle |

Phase-0 recapture (2026-07-15, pre-region tip `2c2470b2`, after the slim
RingGroup landing + cold evictions): RingSpin `rs_ring_scan` -Os 43–57 ms /
-O3 29–39 ms (med 33.3); DisplacementField `df_fused_scan` NOISE dwell -Os
56–62 ms / -O3 41–46 ms — the acceptance comparisons in
`docs/profiles/shipping/` use these same-tip figures, not the table above.
| shared leaf | `filter_blend` per blend (IslamicStars / RingSpin fused) | 192 / 162 cyc | 86 / 107 cyc | — |

RingSpin's -Os figure is per-blend-chain bound (~960 scan cyc/blend, 4×
overdraw; see the 2026-07-15 report) — the fused-tip profile is the direct
measurement that only the optimization level, not further structure work,
closes its gap. If Phase 0 step 4's recapture shows the -Os column has
drifted, update this table in the same commit.

## 8. Verification matrix

| Check | Command / method | Pass criterion |
|---|---|---|
| Gate self-tests | `just teensy-gate-test` | green, incl. new component cases |
| Device size gates | `just teensy-size` | all four envs pass; phantasm RAM1 `free` unchanged from baseline; component ceiling green |
| Holosphere untouched | compare `teensy_size` output (or hex hash) vs baseline | identical |
| Host suite | `just test` | green; goldens unchanged (macros expand empty under clang) |
| WASM | `just smoke` | green |
| restated FP flags reach codegen | three-way disassembly of one region function (`arm-none-eabi-objdump -d`, symbol via the map): (a) vs a throwaway build with bare `optimize("O3")` — identical output means the fast-math umbrella was ignored (§4.2 #1); (b) vs the same symbol in the `profile_o3` image — a match (modulo addresses/layout) proves the region reproduces true global-O3+fast-math codegen. `vfma` presence alone proves **nothing**: GNU mode contracts FMAs by default even at plain -Os; (c) the `filter.h` sink's `isfinite` NaN-guard branch still present (not folded away) | (a) differs, (b) matches, (c) present |
| Pragma actually applied | per-region map diff (§7.2) | expected symbols grow |
| Speed | `just profile <sentinels>` (§7.3) | tier crossings per the kept-region ledger; hot scopes fully covered by regions land within ~15% of the -O3 ceiling figure (§7 sentinel-baselines table). Judge each region by its **own hot-scope counter** (the baselines table lists counters, not frame time): a whole-frame gap while a sibling region is unlanded is expected stacking, not a wall. A per-counter gap >15% after all serving regions land means an -Os wall is still inside the loop — find it before accepting |
| On-device visual sanity | flash phantasm, observe the roster | no artifacts (manual; float codegen inside regions changes → device output may differ at the LSB from the pure -Os image, which is within the sim/device parity contract — that contract is the integer wrap guards, not float bits) |

## 9. Risks and fallbacks

- **Per-effect instantiation blowup (mainly R3/R5).** Detected by the map diff.
  Fallback: narrow the region to the strategy/leaf loops (`HS_O3_FN`), or offset
  by moving more measured-cold code to flash with the existing `HS_COLD` /
  `HS_COLD_MEMBER` recipes (see the phantasm ITCM ledger) — flash is free.
- **Pragma silently not applied** (GCC template/lambda edge cases). Detected by
  a ~0 size **and** ~0 profile delta; fallback `HS_O3_FN` per function.
- **An -Os wall inside a loop** (helper not covered, not `always_inline`).
  Detected by the ≥15% gap rule in §8; fix by extending the region to the
  helper, never by marking hot helpers `always_inline` solely to smuggle them
  in (that changes the -Os build too).
- **A region on the revert bubble.** Before reverting a region whose win is
  real but whose ITCM cost breaks the budget, A/B an -O2 variant of it: -O3's
  marginal size over -O2 (aggressive inlining/unrolling) is exactly the ITCM
  problem, and on the scalar Cortex-M7 (no SIMD) much of that extra work buys
  little. The owner's rejection of *global* -O2 does not cover this per-region
  use. If an -O2 region is kept, add `HS_O2_BEGIN`/`HS_O2_END` twins to the
  §4.1 block (same flag restatement, level `"O2"`) rather than open-coding a
  divergent pragma.
- **Future code growth eats the padding.** That is the point of the component
  ceiling: it turns the invisible drift into a reviewed one-line budget edit,
  exactly like every other ceiling in `teensy_budgets.json`.

## 10. Acceptance checklist

- [ ] §6 gate extension landed first, with tests, doc update, and a biting
      calibrated ceiling.
- [ ] `HS_O3_BEGIN/END/FN` in `platform.h`, matching §4.1 exactly (flag
      restatement included), no-op everywhere but device-GCC `-Os`.
- [ ] Regions from §5 applied in order, one commit each (each kept region's
      commit bumping the component ceiling, §6.4), ledger filled in (including
      reverted regions and why).
- [ ] ITCM bank count unchanged from Phase 0, with ≥ 4,096 B intra-bank
      padding; RAM1 `free` unchanged.
- [ ] Full §8 matrix green; tier crossings demonstrated on device for at least
      MeshFeedback, Flyby, and Voronoi (the three full 8→16 candidates), or the
      ledger documents the measured reason a crossing was not achievable within
      budget.
