# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/voronoi_ship.log`, captured 2026-10-23 13:33 America/Los_Angeles on COM4.
Replaces `profile_voronoi_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Voronoi, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Voronoi profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:50216, data:147068, headers:8540   free for files:1825792` / `RAM1: variables:312960, code:10344, padding:22424   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `a2e1ea4c43ba0ebe16fdc61146101c3f4428b7ea58e11cf04f51f34e73096090`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **9.242 ms** (frame 2), mean **8.425 ms**, **0/447 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–448 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **13.222 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 457 raw frame rows, including that trailing portion, the observed peak is 13.222 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `vo_shade` averages 7.232 ms/frame; its worst window averages 7.854 ms/frame (frames 417–448).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.237 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.148 | 0.97 | 0.25% |
| ↳ `vo_shade` DUPLICATE-NAME | 0.393 | 0.03 | 0.66% |
| ↳ `vo_kdtree` | 0.343 | 1.00 | 0.58% |
| ↳ `vo_animate` | 0.043 | 1.00 | 0.07% |
| ↳ `canvas_clear` | 0.090 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 50.597 | 1.00 | 85.41% |

Wall min/mean/max: 9.243/59.251/62.971 ms.

### Worst dominant-scope window (frames 417–448)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.448 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.152 | 1.00 | 0.24% |
| ↳ `vo_shade` | 7.854 | 1.00 | 12.58% |
| ↳ `vo_kdtree` | 0.345 | 1.00 | 0.55% |
| ↳ `vo_animate` | 0.042 | 1.00 | 0.07% |
| ↳ `canvas_clear` | 0.085 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 53.950 | 1.00 | 86.39% |

Wall min/mean/max: 62.096/62.462/62.820 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2382.821 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1372.59 | 0.000438 | 0.001650 | 0.014336 | 3.04% |
| `isr_pack` | 135.66 | 0.008740 | 0.009216 | 0.011725 | 1.67% |
| `isr_dma_submit` | 135.66 | 0.000616 | 0.000935 | 0.001196 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `vo_shade` — 0.393 ms/frame, 0.66% of this root window.
2. `vo_kdtree` — 0.343 ms/frame, 0.58% of this root window.
3. `pov_preserve_half` — 0.148 ms/frame, 0.25% of this root window.

README cells: peak 🟢 9.242, spilled 🟢 0/447 (0.00%).

README setup frame 1: 13.222 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=Voronoi  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (Voronoi 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): vo_shade
  [PASS] scope 'vo_shade' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (232.2 ppm, frames 417-448)

INVALID: Voronoi
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (232.2 ppm, frames 417-448). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Voronoi`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
