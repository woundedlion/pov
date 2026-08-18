# Device/Host Divergence — Coverage Ledger

*A tracked inventory of every place engine behavior forks on a `CORE_TEENSY`/`ARDUINO`-only
constant, each tagged with whether an automated **device-value** test build reaches it.*

## Why this document exists

The simulator is built to *predict* the hardware, not approximate it: identical RNG seed
(`Pcg32(1337)`), bit-exact integer mocks, fixed-point round-trips. The few places that
genuinely fork on a device-only constant are the engine's blind spot, because a normal host
build compiles the host value and the device path ships with **zero automated coverage on the
one platform (Teensy) that has no debugger and no console.**

Each fork carries its own note at the source or in the test that pins it — `core/render/filter.h`,
`tests/test_solids.h`, `tests/test_filter.h`, and `tests/test_h_offset_renorm.h` all say so locally.
This ledger collects them in one place so the question *"which forks are reached by a device-value
test, and which are not?"* has a single, auditable answer.

**Citation style:** rows name the *symbol* (`H_OFFSET`, `addmod8`, `case_empty_fn_call`) and the
file that holds it, never a line number. Line numbers in a ledger that is only worth its
auditability drift silently — `tools/docs_check.py` validates the file half of a `path:line`
span and cannot check the number. Grep for the symbol.

A **device-value test** is a host executable that compiles or asserts against the *hardware*
value of the forking constant — either by recompiling the engine in its own translation unit
with the device value forced (the `fastmath_clamp_check` / `h_offset_renorm_check` tactic), or
by pinning the device constant directly in an assertion. A green ledger row means such a test
exists; a red row is a real device-only path that no automated test currently reaches.

## Ledger

| # | Fork | Device vs host | Behavioral? | Device-value test reaches it? | Tracked by |
|---|------|----------------|-------------|-------------------------------|------------|
| 1 | **`beat88` phase arithmetic** (`beat88`, `core/engine/platform_arduino_mocks.h`) | `(millis-timebase)*bpm88*280 >> 16` runs in 64-bit on the LP64 native test build; the device (FastLED) and wasm32 wrap the product mod 2³² before the shift. | **No** — the result is narrowed to `uint16_t`, i.e. bits [16,31] of the product. The device's mod-2³² wrap only discards bits ≥ 32, so bits [16,31] are identical to the 64-bit result; the host's larger `millis` enters only via its low 32 bits (`P mod 2³² = (low32(millis)·K) mod 2³²`), exactly what the device sees. | ⚪ **N/A** — no behavioral fork. Verified equal across many `millis` values; a non-vacuity guard confirmed no value makes the two paths differ. | Nothing to track: not a behavioral fork. |
| 2 | **`H_OFFSET` latitude offset** (`H_OFFSET`, `core/engine/platform.h` — 3 device / 0 host, `HS_TEST_H_OFFSET`-overridable) | Device appends 3 virtual sub-pole rows: `H_VIRT = H + H_OFFSET`. Threads through `core/math/geometry.h`, `core/math/spherical_field.h`, `core/render/scan.h`, `core/render/plot/cull.h`, `core/render/plot/raster.h`, `core/render/sdf/common.h`, `core/render/sdf/csg.h`, `core/render/filter_feedback.h`, `Dynamo` and `ShapeShifter`. Marquee fork: `Filter::Screen::AntiAlias::plot`'s south-pole Y-clip bilinear-tap renorm fold (`core/render/filter.h`), the device's most numerically subtle hot-loop path. | **Yes** — clips (does not stretch) the image where the LEDs stop short of the south pole, and renormalizes the bilinear tap. | ✅ **Yes.** `tests/h_offset_renorm_check.cpp` recompiles the whole engine with `-DHS_TEST_H_OFFSET=3` and runs the renorm fold against an energy-conservation oracle; `tests/test_geometry.h` injects `OFF = 3` for the `y_to_phi`/`pixel_to_vector` mapping. | Closed by the device-value `h_offset_renorm` build. |
| 3 | **Conway/SolidBuilder scratch budget** (`IslamicStars::init`, `SPLIT_SCRATCH_A_DEFAULT` / `SPLIT_SCRATCH_B_DEFAULT`, `effects/IslamicStars.h`) | `IslamicStars::init` sizes `scratch_arena_a/b` at a 116 KB / 74 KB split via `configure_arenas`. An over-budget recipe traps as a **device-only OOM**. | **Yes** — an over-budget recipe edit OOM-traps only on the device. | ✅ **Yes.** `tests/test_solids.h`'s `check_high_water_for_recipe` asserts the host high-water mark against the **real 116 KB / 74 KB device budget**, and the 64-bit host figure is a conservative *upper* bound on the device's 32-bit-pointer footprint, so a host pass guarantees device fit. | Closed by the real-budget high-water guard. |
| 4 | **Hardware I/O layer** (`#ifdef ARDUINO` guards in `hardware/dma_led.h`, `hardware/dma_led_controller.h`, `hardware/pov_single.h`, `hardware/pov_segmented.h`; `hardware/hd107s_frame.h` stubs `arm_dcache_flush` under `#ifndef ARDUINO`) | eDMA setup, register pokes, and the DMA-wedge watchdog exist only on the device. | **Yes** for the device-side logic; the raw register I/O is not host-observable. | ✅/⚪ **Logic yes, I/O no.** The sync flywheel, segment/column mapping, and epoch scheduler are ported to host and exercised by the `pov_sync` multi-board simulator and the `pov_single`/`pov_segmented` tests. The bare DMA/register writes are device-only by nature and not host-reachable. | Logic closed; register I/O is structurally untestable on host. |
| 5 | **Legacy FastLED `random8`/`random16`** (`random8` / `random16`, `core/engine/platform_arduino_mocks.h`) | Host routes through `Pcg32`; device uses FastLED's LCG. The two streams **intentionally do not match.** | **Yes**, but **divergent by design** — only legacy effects call these; modern effects use `hs::rand_*`, which *do* mirror the device `Pcg32(1337)`. | ⚪ **N/A.** No determinism contract to test; documented per-platform. | Documented at the call site; no contract to track. |
| 6 | **`HS_OS_CYCLES()`** (`HS_OS_CYCLES`, `core/engine/platform_diagnostics.h` — `ARM_DWT_CYCCNT` under `CORE_TEENSY`, `0` host) | `ARM_DWT_CYCCNT` (device) / `0` (host). | **No** — a profiling timestamp; never affects render output. | ⚪ **N/A.** | Non-behavioral. |
| 7 | **`StaticCircularBuffer` host-only narrowing** (the `head`/`tail`/`count` members, `core/engine/static_circular_buffer.h`) | `head`/`tail`/`count` are `uint32_t` rather than `size_t`, so pooled structs keep one layout on both targets. On the device `size_t` *is* `uint32_t`, making the narrowing host-only; identical codegen and behavior on hardware. | **No.** | ⚪ **N/A.** | Non-behavioral. |
| 8 | **Packed `uqadd16` saturating blend-add** (`inline_uqadd16` and `pixel_blend_add_packed`, `core/color/color.h`) | Device runs the ARM `uqadd16` DSP instruction (`__ARM_FEATURE_DSP` inline asm); the host runs a portable per-lane software model. Both are specified as two independent 16-bit unsigned saturating adds (g\|b packed in one word, r in another). `lerp16` is **not** a fork — it is portable C on both platforms (signed-multiply asm is avoided, see its body in `core/color/color.h`), so it stays bit-identical. | **Yes** — the packed result feeds `Pixel::operator+=` / `blend_add`. | ✅/⚪ **Model yes, instruction no.** `tests/test_color.h`'s `test_blend_add_packed_lane_layout` compiles the software `uqadd16` plus the exact device lane-packing and pin them against an independent per-channel saturating reference, so the layout the asm relies on is covered. The literal `uqadd16` instruction never executes in CI (no device build), but it implements the same fixed ISA semantics as the model, so a divergence would require an asm typo / ISA misread — structurally untestable on host, like row 4's register I/O. | This row. |
| 9 | **Empty `Fn` invoke** (the `Fn` alias, `core/engine/platform.h`) | Host/WASM `Fn` is `hs::inplace_function`, whose empty-state invoke fail-fast traps through `check_fail`. The device `Fn` is the vendored `teensy::inplace_function`, whose Teensy-core header defines `SG14_INPLACE_FUNCTION_THROW(x)` as `return static_cast<R>(0)` before its own `#ifndef` guard, so the app layer cannot override it: an unbound call returns a zero-initialized `R` and execution continues. | **Yes**, on the bug path only — an unbound `Fn` is a programmer error either way, but the host aborts loudly while the one target with no console silently yields zero. Bound calls are identical. | ❌ **No.** `tests/test_death.h`'s `case_empty_fn_call` proves the trap for `hs::inplace_function`; the death harness is host-only, so nothing exercises the device backend. Closing it means routing the device through `hs::inplace_function` too, which changes device codegen and ITCM footprint (`docs/ledgers/itcm_ledger.md`) for a bug-path-only gain. | This row. |
| 10 | **Pole-LOD knob mutability** (`pole_lod_aggressiveness` / `POLE_LOD_ENABLED`, `core/engine/constants.h`) | Under `ARDUINO`, `pole_lod_aggressiveness` is a `constexpr float` fixed at `HS_POLE_LOD_DEFAULT` (0.0f — no build env overrides it) and `POLE_LOD_ENABLED` folds to `HS_POLE_LOD_DEFAULT > 0.0f`, i.e. **false** on the shipping firmware, so every `if constexpr (POLE_LOD_ENABLED)` guard in `core/render/scan.h` compiles out and firmware has no setter. Host and WASM declare the knob as a mutable global and pin `POLE_LOD_ENABLED` true, so the decimation is always compiled in and gated at runtime by `Scan::pole_lod_run` (settable live via `setPoleLod`). | **No** at equal knob values — `pole_lod_run` returns 1 at aggressiveness 0, so every block offer is one column and the walk is bit-identical to the device's compiled-out walk. The fork is codegen and settability (the device value is frozen at build time), not image content. | ⚪ **N/A / host-unrepresentable.** No host build can compile the `POLE_LOD_ENABLED == false` state, so the device's *shape* is out of reach of any host test; its *behavior* is not, and `tests/test_scan.h` pins both ends of the knob (0 bit-identical to an undecimated walk; 1.0 canvas-anchored runs). A firmware build with `HS_POLE_LOD_DEFAULT > 0` would run exactly the code those tests drive, with the value folded in. | This row. Revisit if firmware gains a runtime setter or ships a non-zero `HS_POLE_LOD_DEFAULT`. |
| 11 | **`addmod8` zero modulus** (`addmod8`, `core/engine/platform_arduino_mocks.h`) | FastLED's non-AVR `addmod8` is `a += b; while (a >= m) a -= m;` (`lib8tion/math8.h`), so `m == 0` subtracts zero forever and hangs the device. The host mock returns the `uint8_t`-wrapped sum instead. Every `m != 0` input is bit-identical, including the past-255 wrap (`addmod8(200, 100, 7)` is `44 % 7`, not `300 % 7`). | **Yes**, on the bug path only — a zero modulus is a caller error either way, but the one target with no console hangs while the host returns and continues. No shipped call site reaches it: every caller in `core/engine/effects_legacy.h` passes a nonzero dimension or literal (`W`, `H`, `16`, `sizeof(counts)`). | ❌ **No.** `tests/test_platform.h`'s `test_addmod8_wraps_before_reducing` pins the host side including `m == 0`, and sweeps all 256×256 addends against a repeated-subtraction reference at `m = 7`. A non-terminating device loop is not observable from a host test; closing it means a guard in FastLED's own header. | This row. |

Legend: ✅ reached by a device-value test · ❌ real device-only path with no device-value test
· ⚪ no behavioral fork / divergent by design (nothing to cover).

## Standing risk

**Two red rows, both bug-path-only, both accepted risks.** Row 9 (empty `Fn` invoke) fires only
when a caller invokes an unbound callable; the cost of closing it is device codegen and ITCM
footprint on every `Fn` instantiation. Row 11 (`addmod8` zero modulus) fires only on a zero
modulus, which no shipped call site passes; closing it means guarding inside FastLED's own
header. Both share the shape that makes them tolerable — the device diverges only where the
caller is already wrong — and the shape that makes them worth recording: the device failure is
silent (a hang, a zero return) on the one target with no console.

Row 1 (`beat88`) is not a breach: the `uint16_t` result extracts bits [16,31] of the phase product,
which the device's mod-2³² wrap cannot change, so the 64-bit native build and the 32-bit
device/wasm builds produce bit-identical phases. Every behavioral fork other than the two red
rows above either has device-value coverage (rows 2–4), host model-level coverage with only a
device-only ISA-instruction tail (row 8), or is divergent by design / non-behavioral (rows 5–7,
10).

## Maintenance rule

When a new behavior is gated on a `CORE_TEENSY`/`ARDUINO`-only constant — or an existing
device-only constant grows a new dependent path — add a row here and state plainly whether a
device-value test reaches it. A new fork with no green test is a known, accepted risk only once
it is written down; until then it is the kind of silent sim≠device gap this ledger exists to
surface.
