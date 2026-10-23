# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/hyperlattice_ship.log`, captured 2026-10-23 14:18 America/Los_Angeles on COM4.
Replaces `profile_hyperlattice_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 300000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh HyperLattice profile 300 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:144464, data:152396, headers:8292   free for files:1726464` / `RAM1: variables:312992, code:12888, padding:19880   free for local variables:178528` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `49d6c54a5e79fd24fc8158f2caee215f81076c6a177685f1be3598941c4eaf3e`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **55.035 ms** (frame 3209), mean **30.520 ms**, **0/4768 spilled frames** (0.00%).
The frame scope is 1–4768; no warmup or setup frame is excluded.
Frame 1 renders 43.131 ms. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 4776 raw frame rows, including that trailing portion, the observed peak is 55.035 ms and 0 frames spill.

`hl_shader_draw` averages 27.934 ms/frame; its worst window averages 44.074 ms/frame (frames 3089–3104).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 58.901 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 0.94 | 0.24% |
| ↳ `hl_shader_draw` DUPLICATE-NAME | 2.656 | 0.06 | 4.51% |
| ↳ `hl_timeline_step` | 0.008 | 1.00 | 0.01% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 32.226 | 1.00 | 54.71% |

Wall min/mean/max: 24.872/58.904/64.198 ms.

### Peak-render window (frames 3201–3216)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.048 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 1.00 | 0.22% |
| ↳ `hl_shader_draw` | 43.472 | 1.00 | 70.06% |
| ↳ `hl_timeline_step` | 0.040 | 1.00 | 0.06% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 16.116 | 1.00 | 25.97% |

Wall min/mean/max: 51.524/62.053/74.358 ms.

### Worst dominant-scope window (frames 3089–3104)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.466 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 1.00 | 0.22% |
| ↳ `hl_shader_draw` | 44.074 | 1.00 | 70.56% |
| ↳ `hl_timeline_step` | 0.021 | 1.00 | 0.03% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 15.956 | 1.00 | 25.54% |

Wall min/mean/max: 59.762/62.470/65.010 ms.

### Ownership-change window (frames 2673–2688)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 63.196 | 1.00 | 100.00% |
| ↳ `hl_shader_draw` DUPLICATE-NAME | 37.834 | 0.94 | 59.87% |
| ↳ `pov_preserve_half` | 0.139 | 1.00 | 0.22% |
| ↳ `hl_timeline_step` | 0.038 | 1.00 | 0.06% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 20.266 | 1.00 | 32.07% |

Wall min/mean/max: 59.505/63.200/75.338 ms.

### Per-preset table

Ownership buckets include transitions and setup. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `hl_shader_draw` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 6 | — | 55.035 | 0/559 | 44.074 | 34/35 |
| 1 | — | 43.131 | 0/319 | 23.862 | 18/20 |
| 3 | — | 39.202 | 0/559 | 29.124 | 34/35 |
| 9 | — | 37.169 | 0/416 | 33.639 | 25/26 |
| 8 | — | 34.334 | 0/679 | 28.884 | 42/42 |
| 4 | — | 33.722 | 0/439 | 28.590 | 26/27 |
| 5 | — | 32.951 | 0/679 | 29.976 | 42/43 |
| 2 | — | 30.483 | 0/679 | 25.773 | 42/43 |
| 7 | — | 30.042 | 0/439 | 24.817 | 26/27 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 994.348 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1145.56 | 0.000620 | 0.001786 | 0.013870 | 3.29% |
| `isr_pack` | 143.19 | 0.005990 | 0.007748 | 0.011293 | 1.78% |
| `isr_dma_submit` | 143.19 | 0.000770 | 0.000931 | 0.001073 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `hl_shader_draw` — 43.472 ms/frame, 70.06% of this root window.
2. `pov_preserve_half` — 0.139 ms/frame, 0.22% of this root window.
3. `canvas_clear` — 0.089 ms/frame, 0.14% of this root window.

README cells: peak 🟢 55.035 (9), spilled 🟢 0/4768 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=HyperLattice  windows=298
  [PASS] captured >=3 windows (298)
  [PASS] every window names one effect and resolution (HyperLattice 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 298 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 298 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): hl_shader_draw
  [PASS] scope 'hl_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (59.0 ppm, frames 3089-3104)

INVALID: HyperLattice
```

Validation limitations: cycle returns to its first index (peak idx 9, 8 distinct); all 9 presets visited (8); root cyc/600 vs wall sum within 5 ppm (59.0 ppm, frames 3089-3104). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
