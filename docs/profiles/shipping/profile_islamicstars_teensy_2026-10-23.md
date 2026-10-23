# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-10-23, **selective -O3**)

Point-in-time snapshot; reproduce with the command below.
Raw capture: `build/prof/islamicstars_ship.log`, captured 2026-10-23 13:57 America/Los_Angeles on COM4.
Replaces `profile_islamicstars_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions; default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars, 288×144, single-entry playlist, source `49bb199f40a655c23983e73ef6739d47ceb7908c` |
| Method | 306000 ms capture, 16-frame counter windows; shipping show-segment duration |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh IslamicStars profile 306 16` |

The default profile epoch is 3600000 ms. Capture plus 20000 ms attach slack fits inside that epoch.

Profile image size: `FLASH: code:196776, data:207768, headers:9152   free for files:1617920` / `RAM1: variables:313664, code:31016, padding:1752   free for local variables:177856` / `RAM2: variables:520064  free for malloc/new:4224`.

Full Phantasm image size: `FLASH: code:995584, data:851824, headers:9104   free for files:175104` / `RAM1: variables:312736, code:171624, padding:24984   free for local variables:14944` / `RAM2: variables:520064  free for malloc/new:4224`.

Profile ELF SHA-256: `cab38ced4067a98392cc1867247e871d87401fc7801cd4a4ebc20229842e0740`. Phantasm ELF SHA-256: `6d9d66ed551cd05383430124024816cf21a1cf01d7bf0dc0bbdf4bceaaedc0c0`.
The adjacent provenance, compiler/ABI attestation, envdump, build logs and image artifacts identify this capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **58.966 ms** (frame 2056), mean **24.773 ms**, **0/4863 spilled frames** (0.00%).
The scope is per-frame render telemetry (`frame` minus the measured buffer wait), frames 2–4864 in complete counter windows, with the shipping segment clip. Transitions remain included.
**Setup frame 1:** render **16.644 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from the aggregate and ownership tables.

Across all 4873 raw frame rows, including that trailing portion, the observed peak is 58.966 ms and 0 frames spill.

Across all complete counter windows, including the setup-containing window, `is_timeline_step` averages 24.471 ms/frame; its worst window averages 56.323 ms/frame (frames 2049–2064).

A display window is 62.5 ms. Spilled means render exceeds that budget. Render telemetry subtracts `canvas_buffer_wait`; wall time includes display-sync idle. Measurements use the shipping segment clip.

The capture samples this source’s shipping show segment. Preset coverage and wrap diagnostics below apply only to the observed frames; a show can end before the final departure returns to its first preset.

## Phase-by-phase readout

Scopes are inclusive averages. Nested, MIXED-PARENT and DUPLICATE-NAME rows overlap and are not additive.

### Startup (frames 1–16)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 56.350 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 0.94 | 0.24% |
| ↳ `is_timeline_step` | 10.218 | 1.00 | 18.13% |
| ↳ ↳ `is_draw_shape` | 10.058 | 1.00 | 17.85% |
| ↳ ↳ ↳ `is_mesh_scan` | 9.794 | 1.00 | 17.38% |
| ↳ ↳ ↳ `is_face_offsets` | 0.024 | 1.00 | 0.04% |
| ↳ ↳ ↳ `is_mesh_transform` | 0.018 | 1.00 | 0.03% |
| ↳ ↳ ↳ `mesh_visibility` | 0.041 | 1.00 | 0.07% |
| ↳ `is_ripple_prepare` | 0.003 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.102 | 1.00 | 0.18% |
| ↳ `canvas_buffer_wait` | 45.859 | 1.00 | 81.38% |

Wall min/mean/max: 9.341/56.355/63.770 ms.

### Peak-render window (frames 2049–2064)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `mesh_copy_capture` MIXED-PARENT | 0.358 | 22.19 | 0.58% |
| `mesh_capture` MIXED-PARENT | 2.765 | 5.00 | 4.44% |
| `mesh_classify` MIXED-PARENT | 1.484 | 1.00 | 2.39% |
| `frame` | 62.227 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.138 | 1.00 | 0.22% |
| ↳ `is_timeline_step` | 56.323 | 1.00 | 90.51% |
| ↳ ↳ `is_draw_shape` | 56.124 | 1.00 | 90.19% |
| ↳ ↳ ↳ `is_mesh_scan` | 46.672 | 1.00 | 75.00% |
| ↳ ↳ ↳ `is_face_offsets` | 0.033 | 1.00 | 0.05% |
| ↳ ↳ ↳ `is_mesh_transform` | 2.607 | 1.00 | 4.19% |
| ↳ ↳ ↳ `mesh_visibility` | 1.253 | 1.00 | 2.01% |
| ↳ `is_ripple_prepare` | 0.016 | 1.00 | 0.03% |
| ↳ `canvas_clear` | 0.097 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 5.600 | 1.00 | 9.00% |

Wall min/mean/max: 59.462/62.238/65.338 ms.

### Ownership-change window (frames 2145–2160)

| Scope | ms/frame | Calls/frame | Root share |
|---|--:|--:|--:|
| `frame` | 62.513 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.147 | 1.00 | 0.24% |
| ↳ `is_timeline_step` | 35.227 | 1.00 | 56.35% |
| ↳ ↳ `is_draw_shape` | 35.022 | 1.00 | 56.02% |
| ↳ ↳ ↳ `is_mesh_scan` | 32.304 | 1.00 | 51.68% |
| ↳ ↳ ↳ `is_face_offsets` | 0.592 | 1.00 | 0.95% |
| ↳ ↳ ↳ `is_mesh_transform` | 0.165 | 1.00 | 0.26% |
| ↳ ↳ ↳ `mesh_visibility` | 1.192 | 1.00 | 1.91% |
| ↳ `is_ripple_prepare` | 0.002 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.093 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 27.014 | 1.00 | 43.21% |

Wall min/mean/max: 59.319/62.519/65.593 ms.

### Per-preset table

Ownership frame peaks and spill counts include transitions and exclude setup frame 1. Clean scope-window means retain their stated complete-window scope, including the first window. The clean column uses the parser’s single-owner/modal-call selection; it is a window mean, not a frame peak. Continuous morphs do not necessarily provide an independent static hold. Entry 0 is an unmarked startup owner when the marker scheme cannot infer its name.

| Entry | Geometry / name | Peak ms | Spilled/frames | Clean `is_timeline_step` ms/frame | Clean/all windows |
|---|---|--:|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | V=20 E=30 F=12 I=60 | 58.966 | 0/352 | 56.323 | 21/22 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | V=120 E=180 F=62 I=360 | 51.027 | 0/312 | 44.558 | 18/19 |
| truncatedCube_hk90_meta | V=24 E=36 F=14 I=72 | 50.295 | 0/368 | 42.099 | 22/23 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | V=60 E=90 F=32 I=180 | 49.426 | 0/288 | 46.676 | 17/18 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | V=12 E=30 F=20 I=60 | 41.974 | 0/289 | 39.705 | 18/18 |
| octahedron_hk17_ambo_hk73 | V=6 E=12 F=8 I=24 | 35.848 | 0/280 | 24.814 | 17/18 |
| icosahedron_snub_relax_truncate033_hankin62 | V=12 E=30 F=20 I=60 | 34.829 | 0/288 | 32.754 | 17/18 |
| snubDodecahedron_truncate5d_ambo_dual | V=60 E=150 F=92 I=300 | 34.784 | 0/312 | 32.353 | 18/19 |
| dodecahedron_hk72_ambo_dual_hk20 | V=20 E=30 F=12 I=60 | 34.526 | 0/352 | 31.008 | 21/22 |
| dodecahedron_bevel2_relax_gyro | V=20 E=30 F=12 I=60 | 34.132 | 0/352 | 30.816 | 21/22 |
| truncatedIcosahedron_hk58_chamfer63 | V=60 E=90 F=32 I=180 | 33.479 | 0/248 | 31.868 | 14/15 |
| icosahedron_ambo_truncate033_hankin59 | V=12 E=30 F=20 I=60 | 32.285 | 0/272 | 27.722 | 16/17 |
| rhombicuboctahedron_hk63_ambo_hk63 | V=24 E=48 F=26 I=96 | 31.993 | 0/280 | 30.136 | 17/18 |
| dodecahedron_hk54_ambo_hk72 | V=20 E=30 F=12 I=60 | 27.581 | 0/280 | 26.125 | 16/17 |
| icosidodecahedron_truncate5d_ambo_dual | V=30 E=60 F=32 I=120 | 26.976 | 0/312 | 25.655 | 19/20 |
| dodecahedron_hk62_ambo_hk62 | V=20 E=30 F=12 I=60 | 26.656 | 0/278 | 25.486 | 17/18 |

No cycle wrap is observed. These are show-segment measurements, not a demonstrated complete preset cycle. Initial ownership inferred by the parser can appear in the table even when the strict validator counts only advance markers.

### Per-pixel figures

Blended-pixel counts are not enabled in this timing image. Geometric pixel coverage is not a measured blend count.

## Column-ISR / DMA marshaling cost

Steady-peak counter window; separately timestamped ISR interval 998.766 ms.

| ISR scope | Calls/frame | Min ms | Mean ms | Max ms | CPU share |
|---|--:|--:|--:|--:|--:|
| `isr_wake` | 1150.69 | 0.000485 | 0.002030 | 0.024323 | 3.74% |
| `isr_pack` | 143.81 | 0.008626 | 0.009353 | 0.012015 | 2.15% |
| `isr_dma_submit` | 143.81 | 0.000653 | 0.000935 | 0.009393 | 0.21% |

`isr_wake` includes packing and submission; those shares must not be added. Submit is CPU time and does not measure asynchronous wire-transfer duration. Render scopes already absorb interrupt time. Other interrupts are not separately measured.

## Summary ranking

Hot scopes below describe the complete counter window containing the steady peak frame; they are window means.

1. `is_timeline_step` — 56.323 ms/frame, 90.51% of this root window.
2. `mesh_capture` — 2.765 ms/frame, 4.44% of this root window.
3. `mesh_classify` — 1.484 ms/frame, 2.39% of this root window.

README cells: peak 🟢 58.966 (16), spilled 🟢 0/4863 (0.00%).

README setup frame 1: 16.644 ms. Peak and spilled cells use the steady frame scope.

Source, schedule, duration and exclusions identify this observation. Earlier report peaks are different captures with their own source and frame scope; they are not matched comparisons.

## Caveats

- All cycle scopes include interrupt time and instrumentation overhead.
- Steady figures exclude only frame 1; they include subsequent transitions and are limited to complete counter windows. Counter-scope tables retain their named window scope.
- Global-O3 reference captures retain their own source and scope.

Parser validation:

```text
effect=IslamicStars  windows=304
  [PASS] captured >=3 windows (304)
  [PASS] every window names one effect and resolution (IslamicStars 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] shape markers advance and return to the first (cycle closed): 16 distinct
  [PASS] per-frame render is render, not wall (0 of 304 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 304 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): is_mesh_transform, mesh_capture, mesh_classify, mesh_copy_capture, mesh_sweep_paths, mesh_visibility, scan_face_setup, scan_face_sweep, scan_mesh_raster
  [PASS] scope 'is_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (171.7 ppm, frames 2049-2064)

INVALID: IslamicStars
```

Validation limitations: shape markers advance and return to the first (cycle closed): 16 distinct; root cyc/600 vs wall sum within 5 ppm (171.7 ppm, frames 2049-2064). Direct per-frame render telemetry supports the cadence figures. This report does not claim fully validated cycle-counter exactness or an unobserved wrap.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`. The reproduction command builds, flashes and captures through the board lock. No fixed preset, dwell compression or epoch override is used.
