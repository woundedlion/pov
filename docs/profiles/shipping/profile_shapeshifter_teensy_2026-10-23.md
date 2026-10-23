# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/shapeshifter_ship.log`, captured 2026-10-23 13:38 America/Los_Angeles on COM4.
Replaces `profile_shapeshifter_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 269000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh ShapeShifter profile 269 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:192816, data:155844, headers:8716   free for files:1674240` / `RAM1: variables:312992, code:53832, padding:11704   free for local variables:145760` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `9b442c57bf2f4c86e1ed45e196c0029e90a568374dfaf2505a25c274c3f64403`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **70.224 ms** (frame 3554), mean **39.203 ms**, **134/3999 spilled frames** (3.35%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–4000 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **47.913 ms**, 0/1 spilled. The 12 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 4012 raw frame rows, including that trailing portion, the observed peak is 70.224 ms and 134 frames spill.

Across all complete counter windows, including the setup-containing window, `ss_draw_all` averages 38.820 ms/frame; its worst window averages 66.465 ms/frame (frames 3553–3568).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.660 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.129 | 0.94 | 0.22% |
| ↳ `ss_draw_all` | 34.062 | 1.00 | 57.09% |
| ↳ ↳ `ss_plot_dispatch` | 33.470 | 281.19 | 56.10% |
| ↳ `ss_timeline_step` | 0.071 | 1.00 | 0.12% |
| ↳ `ss_buffer_wait` | 25.383 | 1.00 | 42.55% |
| ↳ ↳ `canvas_clear` | 0.095 | 1.00 | 0.16% |
| ↳ ↳ `canvas_buffer_wait` | 25.281 | 1.00 | 42.37% |

Wall min/mean/max: 34.417/59.674/66.138 ms.

### Peak-render window (frames 3553–3568)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 171.545 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.137 | 1.00 | 0.08% |
| ↳ `ss_draw_all` | 66.465 | 1.00 | 38.75% |
| ↳ `ss_timeline_step` | 0.083 | 1.00 | 0.05% |
| ↳ `ss_buffer_wait` | 104.825 | 1.00 | 61.11% |
| ↳ ↳ `canvas_clear` | 0.091 | 1.00 | 0.05% |
| ↳ ↳ `canvas_buffer_wait` | 104.715 | 1.00 | 61.04% |

Wall min/mean/max: 67.955/171.561/191.932 ms.

### Ownership-change window (frames 3345–3360)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 70.551 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.137 | 1.00 | 0.19% |
| ↳ `ss_draw_all` | 56.313 | 1.00 | 79.82% |
| ↳ `ss_timeline_step` | 0.098 | 1.00 | 0.14% |
| ↳ `ss_buffer_wait` | 13.970 | 1.00 | 19.80% |
| ↳ ↳ `canvas_clear` | 0.098 | 1.00 | 0.14% |
| ↳ ↳ `canvas_buffer_wait` | 13.860 | 1.00 | 19.65% |

Wall min/mean/max: 51.694/70.566/176.526 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `ss_draw_all` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 8 | — | 70.224 | 90/479 | 66.465 | 29/29 |
| 7 | — | 66.797 | 44/479 | 64.123 | 29/30 |
| 1 | — | 57.503 | 0/470 | 47.781 | 29/30 |
| 6 | — | 57.146 | 0/479 | 55.161 | 29/30 |
| 4 | — | 42.311 | 0/479 | 40.782 | 29/30 |
| 5 | — | 37.099 | 0/479 | 31.719 | 29/30 |
| 2 | — | 36.058 | 0/479 | 30.928 | 29/30 |
| 9 | — | 35.983 | 0/176 | 34.397 | 11/11 |
| 3 | — | 29.881 | 0/479 | 25.213 | 29/30 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 2746.957 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 3164.69 | 0.000565 | 0.001933 | 0.023003 | 3.56% |
| `isr_pack` | 395.56 | 0.008053 | 0.008873 | 0.012573 | 2.04% |
| `isr_dma_submit` | 395.56 | 0.000678 | 0.000946 | 0.006466 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `ss_draw_all` — 66.465 ms/frame, 38.75% of this root window.
2. `pov_preserve_half` — 0.137 ms/frame, 0.08% of this root window.
3. `ss_timeline_step` — 0.083 ms/frame, 0.05% of this root window.

README cells: peak 🔴 70.224 (2)<br>🟢 57.503 (7), spilled 🔴 134/958 (13.99%)<br>🟢 0/3041 (0.00%).

README setup frame 1: 47.913 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=ShapeShifter  windows=250
  [PASS] captured >=3 windows (250)
  [PASS] every window names one effect and resolution (ShapeShifter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 250 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 250 windows are short)
  [PASS] scope 'ss_draw_all' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (96.5 ppm, frames 3553-3568)

INVALID: ShapeShifter
```

Validation limitations: cycle returns to its first index (peak idx 9, 8 distinct); all 9 presets visited (8); root cyc/600 vs wall sum within 5 ppm (96.5 ppm, frames 3553-3568). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
