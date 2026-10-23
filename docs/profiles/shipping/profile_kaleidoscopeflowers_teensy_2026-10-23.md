# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/kaleidoscopeflowers_ship.log`, captured 2026-10-23 14:35 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopeflowers_teensy_2026-10-22.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeFlowers, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 173000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeFlowers profile 173 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:88880, data:155232, headers:8816   free for files:1778688` / `RAM1: variables:313024, code:22856, padding:9912   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `8beb455a1028aa779b5d148836432d8a572bf46fb80b530a34457328f9e5b03f`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

Column strobe is off in this shipping configuration.

## Frame cadence

**Complete-window aggregate, including setup:** peak render **53.438 ms** (frame 1), mean **30.585 ms**, **0/2736 spilled frames** (0.00%).
The frame scope is 1–2736; no warmup or setup frame is excluded.
Frame 1 renders 53.438 ms. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 2745 raw frame rows, including that trailing portion, the observed peak is 53.438 ms and 0 frames spill.

`fx_shader_draw` averages 26.905 ms/frame; its worst window averages 29.680 ms/frame (frames 1937–1952).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 59.693 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.132 | 0.94 | 0.22% |
| ↳ `fx_shader_draw` | 25.982 | 1.00 | 43.53% |
| ↳ ↳ `fx_prepare_column` | 0.775 | 144.00 | 1.30% |
| ↳ `fx_prepare_frame` | 1.256 | 1.00 | 2.10% |
| ↳ `fx_advance` | 1.954 | 1.00 | 3.27% |
| ↳ `fx_timeline_step` | 0.110 | 1.00 | 0.18% |
| ↳ `canvas_clear` | 0.104 | 1.00 | 0.17% |
| ↳ `canvas_buffer_wait` | 29.736 | 1.00 | 49.81% |

Wall min/mean/max: 28.597/59.698/62.697 ms.

### Worst dominant-scope window (frames 1937–1952)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.380 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 29.680 | 1.00 | 47.58% |
| ↳ ↳ `fx_prepare_column` | 2.483 | 144.00 | 3.98% |
| ↳ `fx_prepare_frame` | 1.328 | 1.00 | 2.13% |
| ↳ `fx_advance` | 1.905 | 1.00 | 3.05% |
| ↳ `fx_timeline_step` | 0.212 | 1.00 | 0.34% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 28.684 | 1.00 | 45.98% |

Wall min/mean/max: 60.673/62.385/64.086 ms.

### Ownership-change window (frames 593–608)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.438 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 28.514 | 1.00 | 45.67% |
| ↳ ↳ `fx_prepare_column` | 1.303 | 144.00 | 2.09% |
| ↳ `fx_prepare_frame` | 0.628 | 1.00 | 1.01% |
| ↳ `fx_advance` | 1.728 | 1.00 | 2.77% |
| ↳ `fx_timeline_step` | 0.163 | 1.00 | 0.26% |
| ↳ `canvas_clear` | 0.096 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 30.802 | 1.00 | 49.33% |

Wall min/mean/max: 60.554/62.443/64.227 ms.

### Per-preset table

Ownership buckets include transitions and setup. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `fx_shader_draw` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 53.438 | 0/599 | 28.681 | 37/38 |
| 3 | — | 34.712 | 0/1058 | 29.680 | 66/66 |
| 2 | — | 33.518 | 0/1079 | 29.447 | 66/67 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Peak-render counter window; separately timestamped ISR interval 1401.961 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1615.19 | 0.000566 | 0.001523 | 0.030256 | 2.80% |
| `isr_pack` | 130.06 | 0.008716 | 0.009380 | 0.011791 | 1.39% |
| `isr_dma_submit` | 130.06 | 0.000670 | 0.000920 | 0.005426 | 0.13% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

1. `fx_shader_draw` — 25.982 ms/frame, 43.53% of this root window.
2. `fx_advance` — 1.954 ms/frame, 3.27% of this root window.
3. `fx_prepare_frame` — 1.256 ms/frame, 2.10% of this root window.

README cells: peak 🟢 53.438 (3), spilled 🟢 0/2736 (0.00%).

Source, schedule, duration and exclusions identify this observation. It is not a matched comparison with a preceding capture.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- The headline includes setup; no post-warmup peak is presented as a whole-show result.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=KaleidoscopeFlowers  windows=171
  [PASS] captured >=3 windows (171)
  [PASS] every window names one effect and resolution (KaleidoscopeFlowers 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 3, 2 distinct)
  [FAIL] all 3 presets visited (2)
  [PASS] per-frame render is render, not wall (0 of 171 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 171 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (81.9 ppm, frames 1937-1952)

INVALID: KaleidoscopeFlowers
```

Validation limitations: cycle returns to its first index (peak idx 3, 2 distinct); all 3 presets visited (2); root cyc/600 vs wall sum within 5 ppm (81.9 ppm, frames 1937-1952). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeFlowers`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
