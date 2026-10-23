# Comets on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/comets_ship.log`, captured 2026-10-23 13:47 America/Los_Angeles on COM4.
Replaces `profile_comets_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Comets, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 120000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Comets profile 120 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:90152, data:153016, headers:8736   free for files:1779712` / `RAM1: variables:312960, code:26728, padding:6040   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `657cd5f3d60dc4f8d8a19eb308eb474b64aebbc0dab448da0acd236747bd8ed1`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **30.701 ms** (frame 647), mean **15.411 ms**, **0/1888 spilled frames** (0.00%).
The frame scope is 1–1888; no warmup or setup frame is excluded.
Frame 1 renders 0.897 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 1897 raw frame rows, including that trailing portion, the observed peak is 30.701 ms and 0 frames spill.

`cm_draw_trail` averages 14.456 ms/frame; its worst window averages 20.700 ms/frame (frames 513–528).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 54.868 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 0.94 | 0.26% |
| ↳ `cm_draw_trail` | 1.780 | 1.00 | 3.24% |
| ↳ `cm_wipe_rebake` | 0.000 | 1.00 | 0.00% |
| ↳ `cm_timeline_step` | 0.112 | 1.00 | 0.20% |
| ↳ `canvas_clear` | 0.099 | 1.00 | 0.18% |
| ↳ `canvas_buffer_wait` | 52.727 | 1.00 | 96.10% |

Wall min/mean/max: 0.899/54.871/63.457 ms.

### Peak-render window (frames 641–656)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.064 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 1.00 | 0.22% |
| ↳ `cm_draw_trail` | 20.134 | 1.00 | 32.44% |
| ↳ `cm_wipe_rebake` | 2.001 | 1.00 | 3.22% |
| ↳ `cm_timeline_step` | 0.176 | 1.00 | 0.28% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 39.517 | 1.00 | 63.67% |

Wall min/mean/max: 45.519/62.068/77.619 ms.

### Worst dominant-scope window (frames 513–528)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.510 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 1.00 | 0.22% |
| ↳ `cm_draw_trail` | 20.700 | 1.00 | 33.11% |
| ↳ `cm_wipe_rebake` | 2.019 | 1.00 | 3.23% |
| ↳ `cm_timeline_step` | 0.169 | 1.00 | 0.27% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 39.384 | 1.00 | 63.00% |

Wall min/mean/max: 58.207/62.515/66.941 ms.

### Ownership-change window (frames 625–640)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.006 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 1.00 | 0.22% |
| ↳ `cm_draw_trail` | 20.008 | 1.00 | 32.27% |
| ↳ `cm_wipe_rebake` | 0.000 | 1.00 | 0.00% |
| ↳ `cm_timeline_step` | 0.158 | 1.00 | 0.26% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 41.587 | 1.00 | 67.07% |

Wall min/mean/max: 56.984/62.010/66.607 ms.

### Per-preset table

Ownership buckets include transitions and setup. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `cm_draw_trail` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 5 | — | 30.701 | 0/160 | 20.134 | 9/10 |
| 11 | — | 27.641 | 0/160 | 18.162 | 9/10 |
| 10 | — | 27.236 | 0/160 | 20.249 | 9/10 |
| 4 | — | 26.374 | 0/160 | 20.700 | 9/10 |
| 9 | — | 26.063 | 0/160 | 19.134 | 9/10 |
| 2 | — | 25.233 | 0/160 | 19.745 | 9/10 |
| 12 | — | 23.544 | 0/129 | 20.173 | 8/8 |
| 8 | — | 22.527 | 0/160 | 19.538 | 9/10 |
| 6 | — | 21.137 | 0/160 | 15.143 | 9/10 |
| 3 | — | 19.219 | 0/160 | 15.002 | 9/10 |
| 7 | — | 15.905 | 0/160 | 12.607 | 9/10 |
| 1 | — | 14.911 | 0/159 | 12.287 | 9/10 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 994.335 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1145.56 | 0.000610 | 0.001708 | 0.013081 | 3.14% |
| `isr_pack` | 143.19 | 0.005990 | 0.007091 | 0.010680 | 1.63% |
| `isr_dma_submit` | 143.19 | 0.000620 | 0.000936 | 0.001155 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `cm_draw_trail` — 20.134 ms/frame, 32.44% of this root window.
2. `cm_wipe_rebake` — 2.001 ms/frame, 3.22% of this root window.
3. `cm_timeline_step` — 0.176 ms/frame, 0.28% of this root window.

README cells: peak 🟢 30.701 (12), spilled 🟢 0/1888 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=Comets  windows=118
  [PASS] captured >=3 windows (118)
  [PASS] every window names one effect and resolution (Comets 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 12, 11 distinct)
  [FAIL] all 12 presets visited (11)
  [PASS] per-frame render is render, not wall (0 of 118 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 118 windows are short)
  [PASS] scope 'cm_draw_trail' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (86.8 ppm, frames 513-528)

INVALID: Comets
```

Validation limitations: cycle returns to its first index (peak idx 12, 11 distinct); all 12 presets visited (11); root cyc/600 vs wall sum within 5 ppm (86.8 ppm, frames 513-528). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Comets`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
