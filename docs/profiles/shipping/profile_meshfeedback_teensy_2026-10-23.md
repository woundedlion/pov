# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/meshfeedback_ship.log`, captured 2026-10-23 13:51 America/Los_Angeles on COM4.
Replaces `profile_meshfeedback_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MeshFeedback, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 181000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MeshFeedback profile 181 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:192248, data:196628, headers:8436   free for files:1634304` / `RAM1: variables:313120, code:42552, padding:22984   free for local variables:145632` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `cc6fe004e8f752d65a102b93d924b3350def8e1d3f448d9160b054387a7e42cc`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **52.652 ms** (frame 2651), mean **43.662 ms**, **0/2863 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–2864 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **5.659 ms**, 0/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 2872 raw frame rows, including that trailing portion, the observed peak is 52.652 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `mf_feedback_flush` averages 36.561 ms/frame; its worst window averages 40.474 ms/frame (frames 513–528).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 57.553 | 1.00 | 100.00% |
| ↳ `mf_mesh_draw` | 6.096 | 1.00 | 10.59% |
| ↳ `mf_feedback_flush` | 30.922 | 1.00 | 53.73% |
| ↳ ↳ `feedback_composite` | 23.231 | 0.94 | 40.36% |
| ↳ ↳ `feedback_populate` | 7.510 | 0.94 | 13.05% |
| ↳ ↳ `feedback_litscan` | 0.071 | 1.00 | 0.12% |
| ↳ `mf_timeline_step` | 0.126 | 1.00 | 0.22% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `mf_apply_params` | 0.001 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.395 | 1.00 | 0.69% |
| ↳ `canvas_buffer_wait` | 19.967 | 1.00 | 34.69% |

Wall min/mean/max: 5.659/57.559/65.115 ms.

### Peak-render window (frames 2641–2656)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.615 | 1.00 | 100.00% |
| ↳ `mf_mesh_draw` | 9.052 | 1.00 | 14.46% |
| ↳ `mf_feedback_flush` | 32.611 | 1.00 | 52.08% |
| ↳ ↳ `feedback_composite` | 32.016 | 1.00 | 51.13% |
| ↳ ↳ `feedback_populate` | 0.496 | 1.00 | 0.79% |
| ↳ ↳ `feedback_litscan` | 0.003 | 1.00 | 0.00% |
| ↳ `mf_timeline_step` | 0.114 | 1.00 | 0.18% |
| ↳ `mf_apply_params` | 0.043 | 1.00 | 0.07% |
| ↳ `canvas_clear` | 0.394 | 1.00 | 0.63% |
| ↳ `canvas_buffer_wait` | 20.354 | 1.00 | 32.51% |

Wall min/mean/max: 53.731/62.621/74.190 ms.

### Worst dominant-scope window (frames 513–528)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.384 | 1.00 | 100.00% |
| ↳ `mf_mesh_draw` | 7.068 | 1.00 | 11.33% |
| ↳ `mf_feedback_flush` | 40.474 | 1.00 | 64.88% |
| ↳ ↳ `feedback_composite` | 32.317 | 1.00 | 51.80% |
| ↳ ↳ `feedback_populate` | 8.067 | 1.00 | 12.93% |
| ↳ ↳ `feedback_litscan` | 0.003 | 1.00 | 0.00% |
| ↳ `mf_timeline_step` | 0.120 | 1.00 | 0.19% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `mf_apply_params` | 0.001 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.398 | 1.00 | 0.64% |
| ↳ `canvas_buffer_wait` | 14.280 | 1.00 | 22.89% |

Wall min/mean/max: 62.179/62.390/62.602 ms.

### Ownership-change window (frames 1201–1216)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.400 | 1.00 | 100.00% |
| ↳ `mf_mesh_draw` | 7.674 | 1.00 | 12.30% |
| ↳ `mf_feedback_flush` | 39.840 | 1.00 | 63.85% |
| ↳ ↳ `feedback_composite` | 31.720 | 1.00 | 50.83% |
| ↳ ↳ `feedback_populate` | 8.013 | 1.00 | 12.84% |
| ↳ ↳ `feedback_litscan` | 0.003 | 1.00 | 0.01% |
| ↳ `mf_timeline_step` | 0.119 | 1.00 | 0.19% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `mf_apply_params` | 0.017 | 1.00 | 0.03% |
| ↳ `canvas_clear` | 0.398 | 1.00 | 0.64% |
| ↳ `canvas_buffer_wait` | 14.306 | 1.00 | 22.93% |

Wall min/mean/max: 60.525/62.406/64.989 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `mf_feedback_flush` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 12 | — | 52.652 | 0/214 | 33.043 | 13/13 |
| 3 | — | 50.497 | 0/241 | 40.474 | 14/15 |
| 6 | — | 50.388 | 0/241 | 39.278 | 14/15 |
| 11 | — | 49.588 | 0/241 | 32.827 | 14/15 |
| 9 | — | 49.549 | 0/241 | 40.397 | 14/15 |
| 8 | — | 49.545 | 0/241 | 38.909 | 14/15 |
| 5 | — | 49.363 | 0/241 | 40.405 | 14/15 |
| 10 | — | 48.448 | 0/241 | 33.711 | 14/15 |
| 2 | — | 45.197 | 0/241 | 39.355 | 15/16 |
| 7 | — | 44.899 | 0/241 | 36.491 | 14/15 |
| 4 | — | 43.759 | 0/241 | 39.735 | 14/15 |
| 1 | — | 43.686 | 0/239 | 35.179 | 15/15 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 1003.188 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1155.75 | 0.000473 | 0.002058 | 0.023508 | 3.79% |
| `isr_pack` | 144.44 | 0.008338 | 0.009710 | 0.011975 | 2.23% |
| `isr_dma_submit` | 144.44 | 0.000785 | 0.000943 | 0.001810 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `mf_feedback_flush` — 32.611 ms/frame, 52.08% of this root window.
2. `mf_mesh_draw` — 9.052 ms/frame, 14.46% of this root window.
3. `canvas_clear` — 0.394 ms/frame, 0.63% of this root window.

README cells: peak 🟢 52.652 (12), spilled 🟢 0/2863 (0.00%).

README setup frame 1: 5.659 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=MeshFeedback  windows=179
  [PASS] captured >=3 windows (179)
  [PASS] every window names one effect and resolution (MeshFeedback 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 12, 11 distinct)
  [FAIL] all 12 presets visited (11)
  [PASS] per-frame render is render, not wall (0 of 179 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 179 windows are short)
  [PASS] scope 'mf_feedback_flush' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (112.8 ppm, frames 513-528)

INVALID: MeshFeedback
```

Validation limitations: cycle returns to its first index (peak idx 12, 11 distinct); all 12 presets visited (11); root cyc/600 vs wall sum within 5 ppm (112.8 ppm, frames 513-528). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MeshFeedback`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
