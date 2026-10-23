# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/mindsplatter_ship.log`, captured 2026-10-23 13:40 America/Los_Angeles on COM4.
Replaces `profile_mindsplatter_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MindSplatter, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 80000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MindSplatter profile 80 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:152064, data:550400, headers:8192   free for files:1320960` / `RAM1: variables:313152, code:48472, padding:17064   free for local variables:145600` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `37c17839feb4814922b497d80a708327f59ec4732cd4697a5ac45b313a6aa5d1`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **59.437 ms** (frame 1003), mean **27.827 ms**, **0/1247 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–1248 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **0.379 ms**, 0/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 1256 raw frame rows, including that trailing portion, the observed peak is 59.437 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `msp_draw_particles` averages 24.877 ms/frame; its worst window averages 51.028 ms/frame (frames 993–1008).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 54.766 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.140 | 0.94 | 0.26% |
| ↳ `msp_draw_particles` | 0.863 | 1.00 | 1.57% |
| ↳ ↳ `msp_particle_scan` | 0.000 | 0.06 | 0.00% |
| ↳ `msp_particle_step` | 0.135 | 1.00 | 0.25% |
| ↳ `msp_timeline_step` | 0.090 | 1.00 | 0.16% |
| ↳ `canvas_clear` | 0.099 | 1.00 | 0.18% |
| ↳ `canvas_buffer_wait` | 53.395 | 1.00 | 97.50% |

Wall min/mean/max: 0.379/54.770/62.635 ms.

### Peak-render window (frames 993–1008)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.300 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.134 | 1.00 | 0.21% |
| ↳ `msp_draw_particles` | 51.028 | 1.00 | 81.91% |
| ↳ `msp_particle_step` | 2.973 | 1.00 | 4.77% |
| ↳ `msp_timeline_step` | 0.131 | 1.00 | 0.21% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 7.887 | 1.00 | 12.66% |

Wall min/mean/max: 57.992/62.305/66.241 ms.

### Ownership-change window (frames 945–960)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.016 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `msp_draw_particles` | 36.998 | 1.00 | 59.66% |
| ↳ `msp_particle_step` | 7.781 | 1.00 | 12.55% |
| ↳ `msp_timeline_step` | 0.110 | 1.00 | 0.18% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 16.830 | 1.00 | 27.14% |

Wall min/mean/max: 54.175/62.022/68.345 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `msp_draw_particles` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 7 | — | 59.437 | 0/159 | 51.028 | 9/10 |
| 6 | — | 53.666 | 0/159 | 39.469 | 9/10 |
| 2 | — | 50.723 | 0/159 | 41.300 | 9/10 |
| 4 | — | 49.823 | 0/159 | 33.391 | 9/10 |
| 5 | — | 37.570 | 0/159 | 27.009 | 9/10 |
| 8 | — | 33.437 | 0/135 | 28.168 | 8/8 |
| 3 | — | 33.048 | 0/159 | 29.342 | 9/10 |
| 1 | — | 25.679 | 0/158 | 19.493 | 9/10 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 998.828 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1150.75 | 0.000615 | 0.001961 | 0.026830 | 3.61% |
| `isr_pack` | 143.81 | 0.006443 | 0.008490 | 0.011385 | 1.95% |
| `isr_dma_submit` | 143.81 | 0.000770 | 0.000948 | 0.004861 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `msp_draw_particles` — 51.028 ms/frame, 81.91% of this root window.
2. `msp_particle_step` — 2.973 ms/frame, 4.77% of this root window.
3. `pov_preserve_half` — 0.134 ms/frame, 0.21% of this root window.

README cells: peak 🟢 59.437 (8), spilled 🟢 0/1247 (0.00%).

README setup frame 1: 0.379 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=MindSplatter  windows=78
  [PASS] captured >=3 windows (78)
  [PASS] every window names one effect and resolution (MindSplatter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 8, 7 distinct)
  [FAIL] all 8 presets visited (7)
  [PASS] per-frame render is render, not wall (0 of 78 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 78 windows are short)
  [PASS] scope 'msp_draw_particles' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (95.2 ppm, frames 993-1008)

INVALID: MindSplatter
```

Validation limitations: cycle returns to its first index (peak idx 8, 7 distinct); all 8 presets visited (7); root cyc/600 vs wall sum within 5 ppm (95.2 ppm, frames 993-1008). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
