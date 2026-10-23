# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/displacementfield_ship.log`, captured 2026-10-23 13:23 America/Los_Angeles on COM4.
Replaces `profile_displacementfield_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DisplacementField, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 121000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh DisplacementField profile 121 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:117992, data:154896, headers:8712   free for files:1750016` / `RAM1: variables:313472, code:25544, padding:7224   free for local variables:178048` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `67a890b814851ca77d4e4ead9d74b371add78e972fe4e883079affe18368f744`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **51.257 ms** (frame 603), mean **41.391 ms**, **0/1888 spilled frames** (0.00%).
The frame scope is 1–1888; no warmup or setup frame is excluded.
Frame 1 renders 18.953 ms. The 25 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 1913 raw frame rows, including that trailing portion, the observed peak is 51.257 ms and 0 frames spill.

`df_timeline_step` averages 41.089 ms/frame; its worst window averages 49.245 ms/frame (frames 321–352).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.331 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.145 | 0.97 | 0.24% |
| ↳ `df_timeline_step` | 34.539 | 1.00 | 57.25% |
| ↳ ↳ `df_draw_rings` | 34.390 | 1.00 | 57.00% |
| ↳ ↳ ↳ `df_fused_scan` DUPLICATE-NAME | 2.479 | 0.22 | 4.11% |
| ↳ ↳ ↳ `df_lut_bake` DUPLICATE-NAME | 0.605 | 9.00 | 1.00% |
| ↳ ↳ ↳ `df_chunk_cull` DUPLICATE-NAME | 0.135 | 9.00 | 0.22% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 25.515 | 1.00 | 42.29% |
| ↳ `df_prepare_fields` | 0.001 | 1.00 | 0.00% |

Wall min/mean/max: 16.113/60.334/86.367 ms.

### Peak-render window (frames 577–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `ring_stack_table` MIXED-PARENT | 1.819 | 2.00 | 2.91% |
| `frame` | 62.489 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.144 | 1.00 | 0.23% |
| ↳ `df_timeline_step` | 48.225 | 1.00 | 77.17% |
| ↳ ↳ `df_draw_rings` | 48.025 | 1.00 | 76.85% |
| ↳ ↳ ↳ `df_fused_scan` | 25.022 | 1.00 | 40.04% |
| ↳ ↳ ↳ `df_lut_bake` | 16.512 | 42.59 | 26.42% |
| ↳ ↳ ↳ `df_chunk_cull` | 1.574 | 43.91 | 2.52% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 13.986 | 1.00 | 22.38% |
| ↳ `df_prepare_fields` | 0.001 | 1.00 | 0.00% |

Wall min/mean/max: 60.484/62.491/64.824 ms.

### Worst dominant-scope window (frames 321–352)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `ring_stack_table` MIXED-PARENT | 2.123 | 2.00 | 3.40% |
| `frame` | 62.388 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.146 | 1.00 | 0.23% |
| ↳ `df_timeline_step` | 49.245 | 1.00 | 78.93% |
| ↳ ↳ `df_draw_rings` | 49.061 | 1.00 | 78.64% |
| ↳ ↳ ↳ `df_fused_scan` | 25.281 | 1.00 | 40.52% |
| ↳ ↳ ↳ `df_lut_bake` | 16.512 | 47.25 | 26.47% |
| ↳ ↳ ↳ `df_chunk_cull` | 1.758 | 47.62 | 2.82% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.096 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 12.860 | 1.00 | 20.61% |
| ↳ `df_prepare_fields` | 0.001 | 1.00 | 0.00% |

Wall min/mean/max: 61.648/62.390/63.178 ms.

### Ownership-change window (frames 897–928)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.829 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.147 | 1.00 | 0.23% |
| ↳ `df_timeline_step` | 19.777 | 1.00 | 31.48% |
| ↳ ↳ `df_draw_rings` | 19.556 | 1.00 | 31.13% |
| ↳ ↳ ↳ `df_fused_scan` DUPLICATE-NAME | 1.251 | 0.12 | 1.99% |
| ↳ ↳ ↳ `df_lut_bake` DUPLICATE-NAME | 0.281 | 3.69 | 0.45% |
| ↳ ↳ ↳ `df_chunk_cull` DUPLICATE-NAME | 0.063 | 4.38 | 0.10% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.094 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 42.744 | 1.00 | 68.03% |
| ↳ `df_prepare_fields` | 0.006 | 1.00 | 0.01% |

Wall min/mean/max: 58.671/62.832/63.880 ms.

### Per-preset table

Ownership buckets include transitions and setup. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `df_timeline_step` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 51.257 | 0/899 | 49.245 | 28/29 |
| 2 | — | 43.077 | 0/989 | 40.926 | 30/30 |

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2003.328 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1153.97 | 0.000606 | 0.001946 | 0.028798 | 3.59% |
| `isr_pack` | 144.25 | 0.007441 | 0.008778 | 0.014313 | 2.02% |
| `isr_dma_submit` | 144.25 | 0.000628 | 0.000940 | 0.006965 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `df_timeline_step` — 48.225 ms/frame, 77.17% of this root window.
2. `ring_stack_table` — 1.819 ms/frame, 2.91% of this root window.
3. `pov_preserve_half` — 0.144 ms/frame, 0.23% of this root window.

README cells: peak 🟢 51.257 (2), spilled 🟢 0/1888 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=DisplacementField  windows=59
  [PASS] captured >=3 windows (59)
  [PASS] every window names one effect and resolution (DisplacementField 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] cycle returns to its first index (peak idx 2, 2 distinct)
  [PASS] all 2 presets visited (2)
  [PASS] per-frame render is render, not wall (0 of 59 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 59 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): df_chunk_cull, df_fused_scan, df_hue_table_prep, df_lut_bake, ring_stack_table
  [PASS] scope 'df_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (38.5 ppm, frames 321-352)

INVALID: DisplacementField
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (38.5 ppm, frames 321-352). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
