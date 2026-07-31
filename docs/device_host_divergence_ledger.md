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

A **device-value test** is a host executable that compiles or asserts against the *hardware*
value of the forking constant — either by recompiling the engine in its own translation unit
with the device value forced (the `fastmath_clamp_check` / `h_offset_renorm_check` tactic), or
by pinning the device constant directly in an assertion. A green ledger row means such a test
exists; a red row is a real device-only path that no automated test currently reaches.

## Ledger

| # | Fork | Device vs host | Behavioral? | Device-value test reaches it? | Tracked by |
|---|------|----------------|-------------|-------------------------------|------------|
| 1 | **`beat88` phase arithmetic** (`core/engine/platform_arduino_mocks.h:489-491`) | `(millis-timebase)*bpm88*280 >> 16` runs in 64-bit on the LP64 native test build; the device (FastLED) and wasm32 wrap the product mod 2³² before the shift. | **No** — the result is narrowed to `uint16_t`, i.e. bits [16,31] of the product. The device's mod-2³² wrap only discards bits ≥ 32, so bits [16,31] are identical to the 64-bit result; the host's larger `millis` enters only via its low 32 bits (`P mod 2³² = (low32(millis)·K) mod 2³²`), exactly what the device sees. | ⚪ **N/A** — no behavioral fork. Verified equal across many `millis` values; a non-vacuity guard confirmed no value makes the two paths differ. | Nothing to track: not a behavioral fork. |
| 2 | **`H_OFFSET` latitude offset** (`core/engine/platform.h:187` = 3 device / `:1018-1022` = 0 host, `HS_TEST_H_OFFSET`-overridable) | Device appends 3 virtual sub-pole rows: `H_VIRT = H + H_OFFSET`. Threads through `geometry.h`, `scan.h`, `plot.h`, `sdf.h`, `Dynamo`. Marquee fork: `Filter::Screen::AntiAlias::plot`'s south-pole Y-clip bilinear-tap renorm fold (`core/render/filter.h:1190-1197`), the device's most numerically subtle hot-loop path. | **Yes** — clips (does not stretch) the image where the LEDs stop short of the south pole, and renormalizes the bilinear tap. | ✅ **Yes.** `tests/h_offset_renorm_check.cpp` recompiles the whole engine with `-DHS_TEST_H_OFFSET=3` and runs the renorm fold against an energy-conservation oracle; `tests/test_geometry.h:199` injects `OFF=3` for the `y_to_phi`/`pixel_to_vector` mapping. | Closed by the device-value `h_offset_renorm` build. |
| 3 | **Conway/SolidBuilder scratch budget** (`effects/IslamicStars.h:60-62`, split constants at `:161-163`) | `IslamicStars::init` sizes `scratch_arena_a/b` at a 116 KB / 74 KB split via `configure_arenas`. An over-budget recipe traps as a **device-only OOM**. | **Yes** — an over-budget recipe edit OOM-traps only on the device. | ✅ **Yes.** `tests/test_solids.h:609-628` asserts the host high-water mark against the **real 116 KB / 74 KB device budget**, and the 64-bit host figure is a conservative *upper* bound on the device's 32-bit-pointer footprint, so a host pass guarantees device fit. | Closed by the real-budget high-water guard. |
| 4 | **Hardware I/O layer** (`hardware/dma_led.h:24`, `hardware/dma_led_controller.h:49`, `hardware/pov_single.h:21`, `hardware/pov_segmented.h:51`, all `#ifdef ARDUINO`; `hardware/hd107s_frame.h:33` stubs `arm_dcache_flush` under `#ifndef ARDUINO`) | eDMA setup, register pokes, and the DMA-wedge watchdog exist only on the device. | **Yes** for the device-side logic; the raw register I/O is not host-observable. | ✅/⚪ **Logic yes, I/O no.** The sync flywheel, segment/column mapping, and epoch scheduler are ported to host and exercised by the `pov_sync` multi-board simulator and the `pov_single`/`pov_segmented` tests. The bare DMA/register writes are device-only by nature and not host-reachable. | Logic closed; register I/O is structurally untestable on host. |
| 5 | **Legacy FastLED `random8`/`random16`** (`core/engine/platform_arduino_mocks.h:370-401`) | Host routes through `Pcg32`; device uses FastLED's LCG. The two streams **intentionally do not match.** | **Yes**, but **divergent by design** — only legacy effects call these; modern effects use `hs::rand_*`, which *do* mirror the device `Pcg32(1337)`. | ⚪ **N/A.** No determinism contract to test; documented per-platform. | Documented at the call site; no contract to track. |
| 6 | **`HS_OS_CYCLES()`** (`core/engine/platform.h:171-175` = `ARM_DWT_CYCCNT` under `CORE_TEENSY` / `:1003` = `0` host) | `ARM_DWT_CYCCNT` (device) / `0` (host). | **No** — a profiling timestamp; never affects render output. | ⚪ **N/A.** | Non-behavioral. |
| 7 | **`StaticCircularBuffer` host-only narrowing** (`core/engine/static_circular_buffer.h:336-349`) | `head`/`tail`/`count` are `uint32_t` rather than `size_t`, so pooled structs keep one layout on both targets. On the device `size_t` *is* `uint32_t`, making the narrowing host-only; identical codegen and behavior on hardware. | **No.** | ⚪ **N/A.** | Non-behavioral. |
| 8 | **Packed `uqadd16` saturating blend-add** (`core/color/color.h:20-41`; `pixel_blend_add_packed` at `:436-444`) | Device runs the ARM `uqadd16` DSP instruction (`__ARM_FEATURE_DSP` inline asm); the host runs a portable per-lane software model. Both are specified as two independent 16-bit unsigned saturating adds (g\|b packed in one word, r in another). `lerp16` is **not** a fork — it is portable C on both platforms (signed-multiply asm is avoided, see `core/color/color.h:151-154`), so it stays bit-identical. | **Yes** — the packed result feeds `Pixel::operator+=` / `blend_add`. | ✅/⚪ **Model yes, instruction no.** `tests/test_color.h:165`'s `test_blend_add_packed_lane_layout` compiles the software `uqadd16` plus the exact device lane-packing and pin them against an independent per-channel saturating reference, so the layout the asm relies on is covered. The literal `uqadd16` instruction never executes in CI (no device build), but it implements the same fixed ISA semantics as the model, so a divergence would require an asm typo / ISA misread — structurally untestable on host, like row 4's register I/O. | This row. |
| 9 | **Empty `Fn` invoke** (`core/engine/platform.h:1385-1403`) | Host/WASM `Fn` is `hs::inplace_function`, whose empty-state invoke fail-fast traps through `check_fail`. The device `Fn` is the vendored `teensy::inplace_function`, whose Teensy-core header defines `SG14_INPLACE_FUNCTION_THROW(x)` as `return static_cast<R>(0)` before its own `#ifndef` guard, so the app layer cannot override it: an unbound call returns a zero-initialized `R` and execution continues. | **Yes**, on the bug path only — an unbound `Fn` is a programmer error either way, but the host aborts loudly while the one target with no console silently yields zero. Bound calls are identical. | ❌ **No.** `tests/test_death.h:1339`'s `empty_fn_call` case proves the trap for `hs::inplace_function`; the death harness is host-only, so nothing exercises the device backend. Closing it means routing the device through `hs::inplace_function` too, which changes device codegen and ITCM footprint (`docs/itcm_ledger.md`) for a bug-path-only gain. | This row. |

Legend: ✅ reached by a device-value test · ❌ real device-only path with no device-value test
· ⚪ no behavioral fork / divergent by design (nothing to cover).

## Standing risk

**One red row: row 9 (empty `Fn` invoke)** — an accepted risk. It fires only when a caller invokes
an unbound callable, i.e. on a path that is already a bug; the cost of closing it is device
codegen and ITCM footprint on every `Fn` instantiation.

Row 1 (`beat88`) is not a breach: the `uint16_t` result extracts bits [16,31] of the phase product,
which the device's mod-2³² wrap cannot change, so the 64-bit native build and the 32-bit
device/wasm builds produce bit-identical phases. Every other
behavioral fork either has device-value coverage (rows 2–4), host model-level coverage with only a
device-only ISA-instruction tail (row 8), or is divergent by design / non-behavioral (rows 5–7).

## Maintenance rule

When a new behavior is gated on a `CORE_TEENSY`/`ARDUINO`-only constant — or an existing
device-only constant grows a new dependent path — add a row here and state plainly whether a
device-value test reaches it. A new fork with no green test is a known, accepted risk only once
it is written down; until then it is the kind of silent sim≠device gap this ledger exists to
surface.
