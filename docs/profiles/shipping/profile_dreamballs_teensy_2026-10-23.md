# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/dreamballs_ship.log`, captured 2026-10-23 14:00 America/Los_Angeles on COM4.
Replaces `profile_dreamballs_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DreamBalls, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 200000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh DreamBalls profile 200 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:165744, data:196708, headers:8236   free for files:1660928` / `RAM1: variables:313120, code:33320, padding:32216   free for local variables:145632` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `8b5830ba18c3f875055042f2ffa90a590d48140acadd597a7a340eeac7da423a`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **54.311 ms** (frame 2773), mean **29.089 ms**, **0/3167 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–3168 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **49.756 ms**, 0/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 3176 raw frame rows, including that trailing portion, the observed peak is 54.311 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `db_timeline_step` averages 28.818 ms/frame; its worst window averages 50.555 ms/frame (frames 2769–2784).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.249 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.137 | 0.94 | 0.23% |
| ↳ `db_timeline_step` | 38.958 | 1.00 | 64.66% |
| ↳ ↳ `db_draw` | 38.900 | 1.00 | 64.57% |
| ↳ ↳ ↳ `db_draw_scene` | 38.874 | 1.00 | 64.52% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 21.014 | 1.00 | 34.88% |

Wall min/mean/max: 37.478/60.251/64.089 ms.

### Peak-render window (frames 2769–2784)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.346 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 1.00 | 0.23% |
| ↳ `db_timeline_step` | 50.555 | 1.00 | 81.09% |
| ↳ ↳ `db_draw` | 50.441 | 1.00 | 80.91% |
| ↳ ↳ ↳ `db_draw_scene` | 50.411 | 1.00 | 80.86% |
| ↳ `canvas_clear` | 0.087 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 11.517 | 1.00 | 18.47% |

Wall min/mean/max: 57.744/62.348/65.664 ms.

### Ownership-change window (frames 2865–2880)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.378 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.144 | 1.00 | 0.24% |
| ↳ `db_timeline_step` | 44.721 | 1.00 | 75.32% |
| ↳ ↳ `db_draw` | 44.618 | 1.00 | 75.14% |
| ↳ ↳ ↳ `db_draw_scene` | 44.596 | 0.94 | 75.11% |
| ↳ `canvas_clear` | 0.087 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 14.382 | 1.00 | 24.22% |

Wall min/mean/max: 15.664/59.380/64.894 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `db_timeline_step` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 9 | — | 54.311 | 0/320 | 50.555 | 19/20 |
| 1 | — | 47.638 | 0/318 | 44.059 | 19/20 |
| 10 | — | 41.483 | 0/289 | 37.296 | 18/18 |
| 8 | — | 39.877 | 0/320 | 36.916 | 19/20 |
| 4 | — | 31.992 | 0/320 | 28.508 | 19/20 |
| 2 | — | 30.736 | 0/320 | 25.714 | 19/20 |
| 3 | — | 30.502 | 0/320 | 27.606 | 19/20 |
| 7 | — | 26.849 | 0/320 | 24.837 | 19/20 |
| 6 | — | 21.088 | 0/320 | 18.714 | 19/20 |
| 5 | — | 20.018 | 0/320 | 17.632 | 19/20 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 999.468 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1151.50 | 0.000523 | 0.001960 | 0.024375 | 3.61% |
| `isr_pack` | 143.94 | 0.006590 | 0.008048 | 0.010580 | 1.85% |
| `isr_dma_submit` | 143.94 | 0.000785 | 0.000943 | 0.001656 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `db_timeline_step` — 50.555 ms/frame, 81.09% of this root window.
2. `pov_preserve_half` — 0.143 ms/frame, 0.23% of this root window.
3. `canvas_clear` — 0.087 ms/frame, 0.14% of this root window.

README cells: peak 🟢 54.311 (10), spilled 🟢 0/3167 (0.00%).

README setup frame 1: 49.756 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=DreamBalls  windows=198
  [PASS] captured >=3 windows (198)
  [PASS] every window names one effect and resolution (DreamBalls 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 10, 9 distinct)
  [FAIL] all 10 presets visited (9)
  [PASS] per-frame render is render, not wall (0 of 198 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 198 windows are short)
  [PASS] scope 'db_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (46.4 ppm, frames 2769-2784)

INVALID: DreamBalls
```

Validation limitations: cycle returns to its first index (peak idx 10, 9 distinct); all 10 presets visited (9); root cyc/600 vs wall sum within 5 ppm (46.4 ppm, frames 2769-2784). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
