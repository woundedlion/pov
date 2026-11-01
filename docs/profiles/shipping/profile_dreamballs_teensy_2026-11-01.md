# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-11-01, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DreamBalls`).
Raw capture: `build/prof/dreamballs_ship.log`, retained from `C:/work/temp/dreamballs-orbit-profile/candidate.log`, captured 2026-11-01 14:12 America/Los_Angeles on COM4. Replaces `profile_dreamballs_teensy_2026-10-23.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz; COM4, serial 6362390; flywheel + DMA ISRs live |
| Image | `profile` env, `-Os` base + newlib-nano; selective `HS_O3` DreamBalls draw, Plot mesh/column walk/point route and raster kernels; deep profiling off |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master, synchronized timing |
| Effect | DreamBalls 288×144, single-entry playlist, source `bdf02d1ef70f51dfc7c6de8f5a625627d33c90f3` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 230 s capture; all ten presets and return to preset 1; no dwell, epoch or preset-order overrides |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh DreamBalls profile 230 16` |

Compiler: `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Baseline source: `de55152d8d044be96dd539b63ac9063c27a98e5f`.
Profile ELF SHA-256: `3f404a4404b4728c06aa88272146ea1f04005f6c493913102b89676c69d6475a`. Phantasm ELF SHA-256: `524d134dacd458b9fbb6f846d4293aeb25883be67eb7759a5cadd444b167d245`.
Adjacent provenance, envdump, compiler/ABI attestations, build logs and retained ELF artifacts identify both images.

Profile image size: `FLASH: code:162584, data:190564, headers:8324, free:1670144` / `RAM1: variables:314144, code:37048, padding:28488, free:144608` / `RAM2: variables:520064, free:4224`.

Profile deltas against baseline: RAM1 code **+0 B**, RAM1 variables **+0 B**, FLASH data **+0 B**; FLASH code +80 B, headers -80 B.

Full Phantasm image size: `FLASH: code:1073656, data:780772, headers:8228, free:168960` / `RAM1: variables:313760, code:196168, padding:440, free:13920` / `RAM2: variables:520064, free:4224`.

Full Phantasm deltas against baseline: RAM1 code **+0 B**, RAM1 variables **+0 B**, FLASH data **+0 B**; FLASH code +72 B, headers -72 B.

Paired image provenance:

| Capture | Source SHA | Profile ELF SHA-256 | Phantasm ELF SHA-256 |
|---|---|---|---|
| Baseline | `de55152d8d044be96dd539b63ac9063c27a98e5f` | `f63ad5e3bccea75fd8cd9e0d6da02c34822e729d3864cdf33ea6cd24b0c596a5` | `b45f1a7f3ab3d63fbacc2fee1708afa61ae7ece5ebfc7ce3f23c05786688140a` |
| Candidate | `bdf02d1ef70f51dfc7c6de8f5a625627d33c90f3` | `3f404a4404b4728c06aa88272146ea1f04005f6c493913102b89676c69d6475a` | `524d134dacd458b9fbb6f846d4293aeb25883be67eb7759a5cadd444b167d245` |

The full Phantasm size/layout gate and cold warning checks for both Phantasm and Profile passed. The warning checks compiled all five project translation units for each image.

Exactness cross-check: window frames 2769-2784 root counter **599,183,230 cycles / 600 MHz = 998638.716667 us**, versus measured wall sum **998,638 us**, a **0.718 ppm** difference. `parse_profile.py ... validate` also checks the complete capture.

## Frame cadence

**Steady render, frame 1 excluded:** peak **51.563 ms** (frame 2752), mean **29.795 ms**, **0/3647 spilled frames** (0.00%). This uses frames 2-3648 in complete counter windows; transitions remain included.
**Setup frame 1:** render **52.168 ms**, 0/1 spilled. The 9 trailing raw frame rows without a complete counter window are excluded from aggregates and ownership tables.

Across complete windows including setup, peak frame render is 52.168 ms and 0/3648 frames spill. `db_timeline_step` averages 29.529 ms/frame; its worst clean-hold window averages 48.462 ms/frame (frames 2769-2784).

A 62.5 ms display window follows from the 480 RPM half-revolution cadence. The shipping segment clip renders a quadrant, approximately 10,368 pixels. All observed preset ownership buckets hold 16 fps when their spill count is zero. `canvas_buffer_wait` is display-sync round-up idle; per-frame render telemetry subtracts it from wall time. Wall samples can jitter about a display window without implying a render spill.

Matched A/B uses the same frames **2-3200** (3,199 frames), excluding setup and including every transition in that range:

| Image | Mean render ms | Peak render ms | Spilled/frames |
|---|--:|--:|--:|
| Baseline | 29.011 | 51.420 | 0/3199 |
| Safe tighter bound | 28.644 | 51.563 | 0/3199 |

Mean render changed by **-0.367 ms/frame** (1.26% lower mean render time). Peak render changed by **+0.143 ms**; the peak is slightly higher in this pass. This is one paired capture, with a baseline repeat, rather than a statistical significance claim. This is an on-device render-time comparison; both captures use the same board, toolchain, shipping flags and deterministic schedule, with the orbit-bound change as the source difference.


A baseline repeat over identical frames 2-3168 gave 28.984071 versus 28.983767 ms/frame (difference -0.000304 ms); peaks were 51.401 and 51.420 ms. The earlier 200 s log is retained separately and is not substituted for the full-cycle baseline.

## Phase-by-phase readout

Phase schedule: each preset owns approximately 320 frames, with sprite fade/adoption at the boundary; the ten-preset cycle returns to preset 1. Startup, single-owner steady drawing and ownership-change windows are shown separately. Scopes are inclusive, so nested rows are not additive. Tree percentages below use the frame root as denominator.

### Startup (frames 1-16)

```text
frame                       60.477 ms 36.286 Mcyc 100.0%
  pov_preserve_half         139.30 us  83.58 kcyc   0.2% x0.94 148.6us/call
  db_timeline_step          39.677 ms 23.806 Mcyc  65.6%
    db_draw                 39.592 ms 23.755 Mcyc  65.5%
      db_draw_scene         39.587 ms 23.752 Mcyc  65.5%
        db_mesh_plot        38.948 ms 23.369 Mcyc  64.4% x18 2163.8us/call
        db_displace         537.58 us 322.55 kcyc   0.9% x18 29.9us/call
  canvas_clear              102.07 us  61.24 kcyc   0.2% x1 102.1us/call
  canvas_buffer_wait        20.526 ms 12.316 Mcyc  33.9% x1 20526.2us/call
```

Wall min/avg/max = 37.657/60.477/66.324 ms. This window contains setup frame 1 (52.168 ms render) and the first steady draws. Startup is reported separately from the steady/A/B frame scope.

### Worst clean-hold render window (frames 2769-2784)

```text
frame                       62.415 ms 37.449 Mcyc 100.0%
  pov_preserve_half         143.80 us  86.28 kcyc   0.2% x1 143.8us/call
  db_timeline_step          48.462 ms 29.077 Mcyc  77.6%
    db_draw                 48.311 ms 28.987 Mcyc  77.4%
      db_draw_scene         48.298 ms 28.979 Mcyc  77.4%
        db_mesh_plot        46.869 ms 28.122 Mcyc  75.1% x6 7811.5us/call
        db_displace          1.150 ms 690.29 kcyc   1.8% x6 191.7us/call
  canvas_clear               94.60 us  56.76 kcyc   0.2% x1 94.6us/call
  canvas_buffer_wait        13.679 ms  8.208 Mcyc  21.9% x1 13679.2us/call
```

Wall min/avg/max = 59.378/62.415/65.365 ms. This is the highest clean-hold timeline window, selected by single-owner/modal-call criteria. Its mesh draw dominates the render cost; the remaining display-window time is primarily buffer wait.

### Ownership-change window (frames 2865-2880)

```text
frame                       59.587 ms 35.752 Mcyc 100.0%
  pov_preserve_half         146.62 us  87.97 kcyc   0.2% x1 146.6us/call
  db_timeline_step          43.722 ms 26.233 Mcyc  73.4%
    db_draw                 43.557 ms 26.134 Mcyc  73.1%
      db_draw_scene         43.551 ms 26.131 Mcyc  73.1%
        db_mesh_plot        42.236 ms 25.341 Mcyc  70.9% x5.62 7508.6us/call
        db_displace          1.073 ms 643.59 kcyc   1.8% x5.62 190.7us/call
  canvas_clear               95.67 us  57.40 kcyc   0.2% x1 95.7us/call
  canvas_buffer_wait        15.593 ms  9.356 Mcyc  26.2% x1 15592.9us/call
```

Wall min/avg/max = 17.410/59.587/66.324 ms. This window straddles a preset marker and therefore mixes the outgoing and adopted presets. It remains in the per-frame cadence and A/B totals, but is excluded from clean-hold comparisons.

### Per-preset table

Rows are ranked by the worst clean-hold `db_mesh_plot` cost. The clean window is the highest `db_timeline_step` reading among the parser's single-owner/modal `db_draw_scene` windows. Render ms is the mean direct render telemetry in that same window; peak/spill columns cover all owned steady frames, including transitions. Copies are truncated as in production. Blended pixels per frame are not measured in this image. FPS is derived from the selected window wall time; sync jitter produces small deviations around 16.

| # | Source solid | Copies; radius; speed | Blended px/f | Mesh plot ms | Render ms | Peak ms | Spilled/frames | fps | Clean/all windows |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| 9 | Disdyakis Triacontahedron | 6; 0.03; 1 | unmeasured | 46.869 | 48.736 | 51.563 | 0/320 | 16.02 | 19/20 |
| 1 | Rhombicuboctahedron | 18; 0.3; 0.4 | unmeasured | 44.679 | 45.776 | 48.027 | 0/638 | 15.99 | 38/40 |
| 10 | Triakis Icosahedron | 6; 0.03; 1 | unmeasured | 35.155 | 36.366 | 39.677 | 0/320 | 16.00 | 19/20 |
| 8 | Triakis Icosahedron | 6; 0.078; 1 | unmeasured | 34.669 | 35.902 | 38.323 | 0/320 | 16.01 | 19/20 |
| 4 | Icosidodecahedron | 10; 0.16; 1 | unmeasured | 27.168 | 28.074 | 31.535 | 0/320 | 16.03 | 19/20 |
| 3 | Truncated Cuboctahedron | 6; 0.16; 1 | unmeasured | 25.074 | 26.151 | 28.573 | 0/320 | 16.00 | 19/20 |
| 2 | Rhombicosidodecahedron | 6; 0.05; 1 | unmeasured | 24.280 | 25.257 | 28.316 | 0/449 | 16.04 | 27/28 |
| 7 | Triakis Icosahedron | 4; 0.131; 1.89 | unmeasured | 22.892 | 23.906 | 25.463 | 0/320 | 15.97 | 19/20 |
| 6 | Truncated Dodecahedron | 4; 0.179; 1.89 | unmeasured | 16.532 | 17.538 | 19.807 | 0/320 | 16.05 | 19/20 |
| 5 | Snub Cube | 4; 0.153; 2.025 | unmeasured | 16.101 | 16.940 | 18.277 | 0/320 | 16.00 | 19/20 |

Observed preset markers: 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10 → 1 → 2. Return to **1/10 is confirmed**; the capture covers a complete cycle. Clean/all window counts are shown per row. Names and render parameters come from the captured source's `DreamBalls::PRESETS`; repeated solids have distinct radius/speed settings.

### Per-pixel figures

Blended pixel counts, cycles per blend, scan cycles per blended pixel and measured coverage are **not measured**: deep/per-pixel profiling is off. The 10,368-pixel quadrant is geometric area, not a blend count. In a separate native workload run over all ten presets' first 32 frames, midpoint evaluations fell from **177,349 to 18,760 (89.42%)**, and CHORD selection rose from 58.94% to 95.66%; 420,587 assertions passed. These native counts explain the intended planner saving but are not device pixel counts.

## Column-ISR / DMA marshaling cost

The ISR interval associated with frames 2769-2784 is separately timestamped: 999.975 ms. Rates below are calls per rendered frame.

```text
isr_wake       1152.00/f  0.553/1.838/15.015 us  cpu 3.38%
  isr_pack        144.00/f  6.591/8.061/10.746 us  cpu 1.85%
  isr_dma_submit  144.00/f  0.701/0.938/3.856 us  cpu 0.21%
```

- Columns list per-call min/avg/max microseconds. Packing averages 8.061 us/call; CPU submission averages 0.938 us/call. Submission kicks asynchronous transfer and does not time the wire.
- Driver-derived worst-case wire duration is **230 us** for the 600-byte image-plus-black-strobe composite: `HD107SFrame<72>` has 300 bytes/frame; `dma::transfer_us` includes 92 functional-clock ticks per byte at 240 MHz for the 24 MHz SPI request. The image-only duration is 115 us. These are calculated transport bounds, not oscilloscope measurements; the nominal column period is 434.028 us.
- Inclusive `isr_wake` CPU share is 3.38%, leaving approximately **60.387 ms** of foreground CPU in each 62.5 ms display window under this ISR mix. Pack and submit run inside wake; their shares are not added. The measured render scopes already include ISR time, so the observed render margin is compared directly with 62.5 ms, without subtracting ISR a second time. No speedup is required for any zero-spill phase.
- DMA-completion and other interrupts are not separately measured.

## Summary ranking

1. `db_mesh_plot` — 75.09% of the root frame, 46.869 ms/frame in the worst clean-hold window.
2. `db_displace` — 1.84% of the root frame, 1.150 ms/frame in the worst clean-hold window.
3. `pov_preserve_half` — 0.23% of the root frame, 0.144 ms/frame in the worst clean-hold window.

The matched device comparison measures the final safe Bhaskara-aware bounds, not the unsafe ideal-circle prototype. Native planner counts show reduced midpoint work; the device mean change above measures its end-to-end benefit in the shipping configuration.

README cells: peak 🟢 51.563 (10), spilled 🟢 0/3647 (0.00%).

README setup frame 1: 52.168 ms. Peak and spilled cells use the steady frame scope.

## Caveats

- CYCCNT free-runs, so all foreground scopes absorb interrupt time.
- `filter_blend` is not enabled in this timing image. When enabled, it parents under its first caller and can be hidden when that parent has zero calls; no blend counts or per-pixel timing are inferred here.
- Deep/per-pixel scopes are off. Enabling per-pixel scopes adds measurement overhead and would require a separate image for attribution.
- This is the shipping selective-O3 configuration: `HS_HOT_FLASH_MEMBER` draw/capture stages, `HS_O3_FN` timed-edge setup and the Plot column-walk/raster regions run at O3 under the Os translation unit. It is not a global-O3 experiment.
- No dwell compression, fixed preset, epoch stretch or preset-order override was used. The current Profile harness's default one-hour epoch stays beyond the 230 s capture, so no epoch rebuild is expected.
- Both captures use isolated source trees. Their committed source SHAs, compiler, image hashes and memory sizes are recorded above; untracked build/profiling outputs are not source changes.
- Cadence aggregates exclude frame 1 and raw trailing rows without complete counter windows. Per-preset peaks/spills include the transition owned by each frame's marker. Matched A/B fixes the frame range to 2-3200.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, `HS_PROFILE_WINDOW=16`. `just profile DreamBalls` routes build, upload and capture through the shared-device wrapper; the exact paired-capture command is in Setup. The bench board was COM4, serial 6362390, and remained locked through each build/upload/capture operation.
