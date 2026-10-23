# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/kaleidoscopehexoil_ship.log`, captured 2026-10-23 14:25 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopehexoil_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeHexOil, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 105000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeHexOil profile 105 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:88160, data:158784, headers:9056   free for files:1775616` / `RAM1: variables:313024, code:22728, padding:10040   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `34dc4fabcf40eaff976ba54f1666e7ef335d08a06a43d8d57dd58778fb182138`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

Column strobe is off in this shipping configuration.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **69.694 ms** (frame 1), mean **40.496 ms**, **1/1648 spilled frames** (0.06%).
The frame scope is 1–1648; no warmup or setup frame is excluded.
Frame 1 renders 69.694 ms. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 1656 raw frame rows, including that trailing portion, the observed peak is 69.694 ms and 1 frames spill.

`fx_shader_draw` averages 36.863 ms/frame; its worst window averages 39.062 ms/frame (frames 817–832).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 61.456 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.129 | 0.94 | 0.21% |
| ↳ `fx_shader_draw` | 38.017 | 1.00 | 61.86% |
| ↳ ↳ `fx_prepare_column` | 0.706 | 144.00 | 1.15% |
| ↳ `fx_prepare_frame` | 0.812 | 1.00 | 1.32% |
| ↳ `fx_advance` | 2.237 | 1.00 | 3.64% |
| ↳ `fx_timeline_step` | 0.118 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.105 | 1.00 | 0.17% |
| ↳ `canvas_buffer_wait` | 19.839 | 1.00 | 32.28% |

Wall min/mean/max: 40.414/61.467/69.697 ms.

### Worst dominant-scope window (frames 817–832)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.393 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 39.062 | 1.00 | 62.61% |
| ↳ ↳ `fx_prepare_column` | 2.337 | 144.00 | 3.75% |
| ↳ `fx_prepare_frame` | 0.928 | 1.00 | 1.49% |
| ↳ `fx_advance` | 2.143 | 1.00 | 3.43% |
| ↳ `fx_timeline_step` | 0.206 | 1.00 | 0.33% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 19.682 | 1.00 | 31.54% |

Wall min/mean/max: 61.184/62.403/63.732 ms.

### Ownership-change window (frames 593–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.459 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 36.224 | 1.00 | 58.00% |
| ↳ ↳ `fx_prepare_column` | 1.531 | 144.00 | 2.45% |
| ↳ `fx_prepare_frame` | 0.972 | 1.00 | 1.56% |
| ↳ `fx_advance` | 2.062 | 1.00 | 3.30% |
| ↳ `fx_timeline_step` | 0.169 | 1.00 | 0.27% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 22.628 | 1.00 | 36.23% |

Wall min/mean/max: 59.629/62.468/65.070 ms.

### Per-preset table

Ownership buckets include transitions and setup. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `fx_shader_draw` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 69.694 | 1/599 | 38.291 | 37/38 |
| 2 | — | 43.705 | 0/1049 | 39.062 | 65/65 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 1413.757 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1628.75 | 0.000471 | 0.001468 | 0.029943 | 2.70% |
| `isr_pack` | 131.81 | 0.008901 | 0.009506 | 0.011988 | 1.41% |
| `isr_dma_submit` | 131.81 | 0.000731 | 0.000915 | 0.001058 | 0.13% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `fx_shader_draw` — 38.017 ms/frame, 61.86% of this root window.
2. `fx_advance` — 2.237 ms/frame, 3.64% of this root window.
3. `fx_prepare_frame` — 0.812 ms/frame, 1.32% of this root window.

README cells: peak 🔴 69.694 (1)<br>🟢 43.705 (1), spilled 🔴 1/599 (0.17%)<br>🟢 0/1049 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=KaleidoscopeHexOil  windows=103
  [PASS] captured >=3 windows (103)
  [PASS] every window names one effect and resolution (KaleidoscopeHexOil 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 103 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 103 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (165.6 ppm, frames 817-832)

INVALID: KaleidoscopeHexOil
```

Validation limitations: cycle returns to its first index (peak idx 2, 1 distinct); all 2 presets visited (1); root cyc/600 vs wall sum within 5 ppm (165.6 ppm, frames 817-832). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexOil`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
