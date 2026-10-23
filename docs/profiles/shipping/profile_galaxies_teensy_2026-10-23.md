# Galaxies on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/galaxies_ship.log`, captured 2026-10-23 13:17 America/Los_Angeles on COM4.
Replaces `profile_galaxies_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Galaxies, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 60000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Galaxies profile 60 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:91232, data:152588, headers:9108   free for files:1778688` / `RAM1: variables:312960, code:23544, padding:9224   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `ab888b2bb174604f674f3e35525d2c2736231f69a5e2709f010f6985144eaec6`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **45.223 ms** (frame 881), mean **24.959 ms**, **0/928 spilled frames** (0.00%).
The frame scope is 1–928; no warmup or setup frame is excluded.
Frame 1 renders 1.494 ms. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 936 raw frame rows, including that trailing portion, the observed peak is 45.223 ms and 0 frames spill.

`gx_particle_step` averages 15.493 ms/frame; its worst window averages 27.086 ms/frame (frames 865–880).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 54.838 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 0.94 | 0.26% |
| ↳ `gx_draw_particles` | 0.842 | 1.00 | 1.54% |
| ↳ `gx_particle_step` | 0.400 | 1.00 | 0.73% |
| ↳ `gx_timeline_step` | 0.085 | 1.00 | 0.16% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.17% |
| ↳ `canvas_buffer_wait` | 53.260 | 1.00 | 97.12% |

Wall min/mean/max: 1.248/54.839/62.773 ms.

### Peak-render window (frames 881–896)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.322 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.145 | 1.00 | 0.23% |
| ↳ `gx_draw_particles` | 16.579 | 1.00 | 26.60% |
| ↳ `gx_particle_step` | 27.070 | 1.00 | 43.44% |
| ↳ `gx_timeline_step` | 0.085 | 1.00 | 0.14% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 18.344 | 1.00 | 29.43% |

Wall min/mean/max: 60.642/62.323/63.900 ms.

### Worst dominant-scope window (frames 865–880)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.425 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.145 | 1.00 | 0.23% |
| ↳ `gx_draw_particles` | 16.958 | 1.00 | 27.16% |
| ↳ `gx_particle_step` | 27.086 | 1.00 | 43.39% |
| ↳ `gx_timeline_step` | 0.077 | 1.00 | 0.12% |
| ↳ `canvas_clear` | 0.088 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 18.061 | 1.00 | 28.93% |

Wall min/mean/max: 61.882/62.426/63.168 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 998.385 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1150.25 | 0.000623 | 0.002460 | 0.024443 | 4.53% |
| `isr_pack` | 143.75 | 0.006145 | 0.007785 | 0.010268 | 1.79% |
| `isr_dma_submit` | 143.75 | 0.000781 | 0.000943 | 0.010771 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `gx_particle_step` — 27.070 ms/frame, 43.44% of this root window.
2. `gx_draw_particles` — 16.579 ms/frame, 26.60% of this root window.
3. `pov_preserve_half` — 0.145 ms/frame, 0.23% of this root window.

README cells: peak 🟢 45.223, spilled 🟢 0/928 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=Galaxies  windows=58
  [PASS] captured >=3 windows (58)
  [PASS] every window names one effect and resolution (Galaxies 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 58 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 58 windows are short)
  [PASS] scope 'gx_particle_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (19.6 ppm, frames 865-880)

INVALID: Galaxies
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (19.6 ppm, frames 865-880). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Galaxies`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
