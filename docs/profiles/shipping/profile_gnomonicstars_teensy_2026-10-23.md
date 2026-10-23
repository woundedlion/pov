# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/gnomonicstars_ship.log`, captured 2026-10-23 13:24 America/Los_Angeles on COM4.
Replaces `profile_gnomonicstars_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GnomonicStars, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh GnomonicStars profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:75488, data:152040, headers:9016   free for files:1795072` / `RAM1: variables:313024, code:20776, padding:11992   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `c9dd57263014e0f49f38f13a33f899ce12afafd32ac97cd766962485867cedc5`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **26.228 ms** (frame 206), mean **12.286 ms**, **0/448 spilled frames** (0.00%).
The frame scope is 1–448; no warmup or setup frame is excluded.
Frame 1 renders 15.374 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 457 raw frame rows, including that trailing portion, the observed peak is 26.228 ms and 0 frames spill.

`gn_draw_stars` averages 11.920 ms/frame; its worst window averages 18.242 ms/frame (frames 193–224).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.312 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.139 | 0.97 | 0.23% |
| ↳ `gn_draw_stars` | 10.252 | 1.00 | 17.28% |
| ↳ ↳ `gn_star_scan` | 8.512 | 251.75 | 14.35% |
| ↳ `gn_spiral_build` | 0.016 | 0.03 | 0.03% |
| ↳ `gn_timeline_step` | 0.098 | 1.00 | 0.16% |
| ↳ ↳ `animation_mobius_step` | 0.035 | 1.00 | 0.06% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 48.703 | 1.00 | 82.11% |

Wall min/mean/max: 14.903/59.333/72.384 ms.

### Peak-render window (frames 193–224)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.449 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `gn_draw_stars` | 18.242 | 1.00 | 29.21% |
| ↳ ↳ `gn_star_scan` | 16.631 | 238.00 | 26.63% |
| ↳ `gn_timeline_step` | 0.100 | 1.00 | 0.16% |
| ↳ ↳ `animation_mobius_step` | 0.034 | 1.00 | 0.05% |
| ↳ `canvas_clear` | 0.084 | 1.00 | 0.13% |
| ↳ `canvas_buffer_wait` | 43.868 | 1.00 | 70.25% |

Wall min/mean/max: 48.197/62.471/75.689 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 2002.633 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1153.59 | 0.000548 | 0.001728 | 0.025023 | 3.18% |
| `isr_pack` | 144.19 | 0.006126 | 0.007346 | 0.009948 | 1.69% |
| `isr_dma_submit` | 144.19 | 0.000618 | 0.000925 | 0.001270 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `gn_draw_stars` — 18.242 ms/frame, 29.21% of this root window.
2. `pov_preserve_half` — 0.142 ms/frame, 0.23% of this root window.
3. `gn_timeline_step` — 0.100 ms/frame, 0.16% of this root window.

README cells: peak 🟢 26.228, spilled 🟢 0/448 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=GnomonicStars  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (GnomonicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'gn_draw_stars' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (358.7 ppm, frames 193-224)

INVALID: GnomonicStars
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (358.7 ppm, frames 193-224). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
