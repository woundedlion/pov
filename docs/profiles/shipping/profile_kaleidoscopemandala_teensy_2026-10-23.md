# KaleidoscopeMandala on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/kaleidoscopemandala_ship.log`, captured 2026-10-23 14:11 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopemandala_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeMandala, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 105000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeMandala profile 105 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:88824, data:155296, headers:8808   free for files:1778688` / `RAM1: variables:313024, code:22680, padding:10088   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `22125079ac50ff04190c0c0b07aab3b11a8b1aa9e28fcf25b0636e418d400a33`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

Column strobe is off in this shipping configuration.

## Frame cadence

**Steady render, frame 1 excluded:** peak **40.428 ms** (frame 1042), mean **33.749 ms**, **0/1631 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–1632 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **55.678 ms**, 0/1 spilled. The 24 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 1656 raw frame rows, including that trailing portion, the observed peak is 55.678 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `fx_shader_draw` averages 28.914 ms/frame; its worst window averages 32.454 ms/frame (frames 1025–1056).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 61.151 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 0.97 | 0.23% |
| ↳ `fx_shader_draw` | 26.707 | 1.00 | 43.67% |
| ↳ ↳ `fx_prepare_column` | 0.960 | 144.00 | 1.57% |
| ↳ `fx_prepare_frame` | 1.028 | 1.00 | 1.68% |
| ↳ `fx_advance` | 2.312 | 1.00 | 3.78% |
| ↳ `fx_timeline_step` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 30.383 | 1.00 | 49.68% |

Wall min/mean/max: 30.021/61.157/63.698 ms.

### Worst dominant-scope window (frames 1025–1056)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.392 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 32.454 | 1.00 | 52.02% |
| ↳ ↳ `fx_prepare_column` | 1.879 | 144.00 | 3.01% |
| ↳ `fx_prepare_frame` | 4.014 | 1.00 | 6.43% |
| ↳ `fx_advance` | 2.203 | 1.00 | 3.53% |
| ↳ `fx_timeline_step` | 0.201 | 1.00 | 0.32% |
| ↳ `canvas_clear` | 0.090 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 22.958 | 1.00 | 36.80% |

Wall min/mean/max: 60.345/62.398/64.419 ms.

### Ownership-change window (frames 577–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.570 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 29.898 | 1.00 | 47.78% |
| ↳ ↳ `fx_prepare_column` | 1.225 | 144.00 | 1.96% |
| ↳ `fx_prepare_frame` | 1.882 | 1.00 | 3.01% |
| ↳ `fx_advance` | 2.159 | 1.00 | 3.45% |
| ↳ `fx_timeline_step` | 0.117 | 1.00 | 0.19% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 27.915 | 1.00 | 44.61% |

Wall min/mean/max: 60.172/62.575/66.584 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `fx_shader_draw` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 40.428 | 0/1033 | 32.454 | 32/32 |
| 1 | — | 34.242 | 0/598 | 29.550 | 18/19 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2403.334 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1384.41 | 0.000521 | 0.001741 | 0.032643 | 3.21% |
| `isr_pack` | 137.16 | 0.008778 | 0.009433 | 0.014353 | 1.72% |
| `isr_dma_submit` | 137.16 | 0.000501 | 0.000918 | 0.005406 | 0.16% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `fx_shader_draw` — 26.707 ms/frame, 43.67% of this root window.
2. `fx_advance` — 2.312 ms/frame, 3.78% of this root window.
3. `fx_prepare_frame` — 1.028 ms/frame, 1.68% of this root window.

README cells: peak 🟢 40.428 (2), spilled 🟢 0/1631 (0.00%).

README setup frame 1: 55.678 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=KaleidoscopeMandala  windows=51
  [PASS] captured >=3 windows (51)
  [PASS] every window names one effect and resolution (KaleidoscopeMandala 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 51 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 51 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (104.8 ppm, frames 1025-1056)

INVALID: KaleidoscopeMandala
```

Validation limitations: cycle returns to its first index (peak idx 2, 1 distinct); all 2 presets visited (1); root cyc/600 vs wall sum within 5 ppm (104.8 ppm, frames 1025-1056). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeMandala`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
