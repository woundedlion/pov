# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/gsreactiondiffusion_ship.log`, captured 2026-10-23 13:25 America/Los_Angeles on COM4.
Replaces `profile_gsreactiondiffusion_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GSReactionDiffusion, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh GSReactionDiffusion profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:79944, data:343828, headers:8356   free for files:1599488` / `RAM1: variables:313248, code:17224, padding:15544   free for local variables:178272` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `17251223775db76d278dd78c05b29e0ee6d480efa7d2c7d4a9105050b0ea1071`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **48.112 ms** (frame 395), mean **40.902 ms**, **0/447 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–448 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **24.883 ms**, 0/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 456 raw frame rows, including that trailing portion, the observed peak is 48.112 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `grd_render` averages 40.577 ms/frame; its worst window averages 46.985 ms/frame (frames 385–416).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 60.168 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 0.97 | 0.23% |
| ↳ `grd_render` | 24.598 | 1.00 | 40.88% |
| ↳ ↳ `grd_shader_draw` | 9.883 | 6.00 | 16.43% |
| ↳ ↳ `grd_cull_flags` | 0.694 | 6.00 | 1.15% |
| ↳ ↳ `grd_publish` | 1.966 | 6.00 | 3.27% |
| ↳ ↳ `grd_physics` | 5.147 | 6.00 | 8.55% |
| ↳ ↳ `grd_pigment` | 3.078 | 1.00 | 5.12% |
| ↳ ↳ `rd_lattice_basis` | 0.139 | 1.00 | 0.23% |
| ↳ ↳ `grd_color_noise` | 2.806 | 1.00 | 4.66% |
| ↳ `rd_timeline_step` | 0.061 | 1.00 | 0.10% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 35.272 | 1.00 | 58.62% |
| `grd_color_palette` | 0.433 | 0.03 | 0.72% |
| `grd_seed_reaction` | 0.448 | 0.03 | 0.74% |

Wall min/mean/max: 21.413/60.170/64.989 ms.

### Peak-render window (frames 385–416)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.392 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.134 | 1.00 | 0.21% |
| ↳ `grd_render` | 46.985 | 1.00 | 75.31% |
| ↳ ↳ `grd_shader_draw` | 28.354 | 6.00 | 45.45% |
| ↳ ↳ `grd_cull_flags` | 0.701 | 6.00 | 1.12% |
| ↳ ↳ `grd_publish` | 4.126 | 6.00 | 6.61% |
| ↳ ↳ `grd_physics` | 5.167 | 6.00 | 8.28% |
| ↳ ↳ `grd_pigment` | 4.241 | 1.00 | 6.80% |
| ↳ ↳ `rd_lattice_basis` | 0.142 | 1.00 | 0.23% |
| ↳ ↳ `grd_color_noise` | 2.837 | 1.00 | 4.55% |
| ↳ `rd_timeline_step` | 0.052 | 1.00 | 0.08% |
| ↳ `canvas_clear` | 0.087 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 15.126 | 1.00 | 24.24% |

Wall min/mean/max: 61.599/62.393/63.107 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 1998.867 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1151.41 | 0.000505 | 0.002036 | 0.030743 | 3.75% |
| `isr_pack` | 143.91 | 0.007763 | 0.009226 | 0.015670 | 2.12% |
| `isr_dma_submit` | 143.91 | 0.000633 | 0.000943 | 0.007800 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `grd_render` — 46.985 ms/frame, 75.31% of this root window.
2. `pov_preserve_half` — 0.134 ms/frame, 0.21% of this root window.
3. `canvas_clear` — 0.087 ms/frame, 0.14% of this root window.

README cells: peak 🟢 48.112, spilled 🟢 0/447 (0.00%).

README setup frame 1: 24.883 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=GSReactionDiffusion  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (GSReactionDiffusion 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'grd_render' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (27.4 ppm, frames 385-416)

INVALID: GSReactionDiffusion
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (27.4 ppm, frames 385-416). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
