# AlienBrain on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/alienbrain_ship.log`, captured 2026-10-23 14:05 America/Los_Angeles on COM4.
Replaces `profile_alienbrain_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AlienBrain, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 240000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AlienBrain profile 240 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:88568, data:155080, headers:8256   free for files:1779712` / `RAM1: variables:313024, code:22888, padding:9880   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `afadb136d6e71c26763671c342d65cb74956f3d440182ee7047ef786637f97b8`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **26.745 ms** (frame 672), mean **25.739 ms**, **0/3807 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–3808 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **48.180 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 3817 raw frame rows, including that trailing portion, the observed peak is 48.180 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `fx_shader_draw` averages 21.927 ms/frame; its worst window averages 23.212 ms/frame (frames 1–16).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.222 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 0.94 | 0.23% |
| ↳ `fx_shader_draw` | 23.212 | 1.00 | 39.19% |
| ↳ ↳ `fx_prepare_column` | 0.568 | 144.00 | 0.96% |
| ↳ `fx_prepare_frame` | 1.000 | 1.00 | 1.69% |
| ↳ `fx_advance` | 2.133 | 1.00 | 3.60% |
| ↳ `fx_timeline_step` | 0.102 | 1.00 | 0.17% |
| ↳ `canvas_clear` | 0.103 | 1.00 | 0.17% |
| ↳ `canvas_buffer_wait` | 32.169 | 1.00 | 54.32% |

Wall min/mean/max: 25.514/59.224/62.606 ms.

### Worst post-setup scope window (frames 3217–3232)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.402 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 22.496 | 1.00 | 36.05% |
| ↳ ↳ `fx_prepare_column` | 1.034 | 144.00 | 1.66% |
| ↳ `fx_prepare_frame` | 0.806 | 1.00 | 1.29% |
| ↳ `fx_advance` | 2.178 | 1.00 | 3.49% |
| ↳ `fx_timeline_step` | 0.206 | 1.00 | 0.33% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 36.150 | 1.00 | 57.93% |

Wall min/mean/max: 62.089/62.403/62.596 ms.

### Ownership-change window (frames 593–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.432 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 22.211 | 1.00 | 35.58% |
| ↳ ↳ `fx_prepare_column` | 0.807 | 144.00 | 1.29% |
| ↳ `fx_prepare_frame` | 0.979 | 1.00 | 1.57% |
| ↳ `fx_advance` | 2.123 | 1.00 | 3.40% |
| ↳ `fx_timeline_step` | 0.154 | 1.00 | 0.25% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 36.380 | 1.00 | 58.27% |

Wall min/mean/max: 59.108/62.435/65.576 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `fx_shader_draw` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 26.745 | 0/1079 | 22.464 | 66/67 |
| 3 | — | 26.708 | 0/1079 | 22.477 | 67/68 |
| 4 | — | 26.463 | 0/1051 | 22.496 | 65/65 |
| 1 | — | 26.175 | 0/598 | 23.212 | 37/38 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 1399.651 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1612.50 | 0.000528 | 0.001493 | 0.030728 | 2.75% |
| `isr_pack` | 129.75 | 0.007680 | 0.008975 | 0.012308 | 1.33% |
| `isr_dma_submit` | 129.75 | 0.000683 | 0.000935 | 0.007085 | 0.13% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `fx_shader_draw` — 23.212 ms/frame, 39.19% of this root window.
2. `fx_advance` — 2.133 ms/frame, 3.60% of this root window.
3. `fx_prepare_frame` — 1.000 ms/frame, 1.69% of this root window.

README cells: peak 🟢 26.745 (4), spilled 🟢 0/3807 (0.00%).

README setup frame 1: 48.180 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=AlienBrain  windows=238
  [PASS] captured >=3 windows (238)
  [PASS] every window names one effect and resolution (AlienBrain 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 4, 3 distinct)
  [FAIL] all 4 presets visited (3)
  [PASS] per-frame render is render, not wall (0 of 238 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 238 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (43.0 ppm, frames 1-16)

INVALID: AlienBrain
```

Validation limitations: cycle returns to its first index (peak idx 4, 3 distinct); all 4 presets visited (3); root cyc/600 vs wall sum within 5 ppm (43.0 ppm, frames 1-16). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AlienBrain`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
