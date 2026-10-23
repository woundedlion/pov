# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/mermaidskin_ship.log`, captured 2026-10-23 14:20 America/Los_Angeles on COM4.
Replaces `profile_mermaidskin_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MermaidSkin, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 30000 ms capture, 32-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MermaidSkin profile 30 32` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:88888, data:154520, headers:8496   free for files:1779712` / `RAM1: variables:313024, code:25256, padding:7512   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `ad04477652f2a29a1ca707f361d9a2c062953b6ef04b4bbc72eb8073acd4b3da`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **42.073 ms** (frame 193), mean **41.633 ms**, **0/447 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–448 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **70.074 ms**, 1/1 spilled. The 8 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 456 raw frame rows, including that trailing portion, the observed peak is 70.074 ms and 1 frames spill.

Across all complete counter windows, including the setup-containing window, `fx_shader_draw` averages 35.001 ms/frame; its worst window averages 35.805 ms/frame (frames 1–32).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–32)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.006 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.136 | 0.97 | 0.22% |
| ↳ `fx_shader_draw` | 35.805 | 1.00 | 57.74% |
| ↳ ↳ `fx_prepare_column` | 1.224 | 144.00 | 1.97% |
| ↳ `fx_prepare_frame` | 4.059 | 1.00 | 6.55% |
| ↳ `fx_advance` | 1.884 | 1.00 | 3.04% |
| ↳ `fx_timeline_step` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 19.598 | 1.00 | 31.61% |

Wall min/mean/max: 41.461/62.011/70.074 ms.

### Worst post-setup scope window (frames 65–96)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.411 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.144 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 35.079 | 1.00 | 56.21% |
| ↳ ↳ `fx_prepare_column` | 1.304 | 144.00 | 2.09% |
| ↳ `fx_prepare_frame` | 4.082 | 1.00 | 6.54% |
| ↳ `fx_advance` | 1.913 | 1.00 | 3.06% |
| ↳ `fx_timeline_step` | 0.101 | 1.00 | 0.16% |
| ↳ `canvas_clear` | 0.094 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 20.680 | 1.00 | 33.13% |

Wall min/mean/max: 62.131/62.415/62.641 ms.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Setup-containing counter window (the all-frame peak), including frame 1; separately timestamped ISR interval 2415.364 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1391.34 | 0.000503 | 0.001725 | 0.021503 | 3.18% |
| `isr_pack` | 138.00 | 0.008846 | 0.009551 | 0.012103 | 1.74% |
| `isr_dma_submit` | 138.00 | 0.000673 | 0.000940 | 0.003666 | 0.17% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the setup-containing counter window, including frame 1; they are complete-window means. The steady frame peak is reported separately above.

1. `fx_shader_draw` — 35.805 ms/frame, 57.74% of this root window.
2. `fx_prepare_frame` — 4.059 ms/frame, 6.55% of this root window.
3. `fx_advance` — 1.884 ms/frame, 3.04% of this root window.

README cells: peak 🟢 42.073, spilled 🟢 0/447 (0.00%).

README setup frame 1: 70.074 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=MermaidSkin  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (MermaidSkin 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (84.4 ppm, frames 1-32)

INVALID: MermaidSkin
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (84.4 ppm, frames 1-32). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MermaidSkin`, `HS_PROFILE_WINDOW=32`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
