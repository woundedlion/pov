# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-09-30, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).

Raw capture: preserved capture (supporting artifact removed), captured 2026-09-30 19:21 on COM4. Replaces the earlier 2026-09-30 report from before pigment blending, hue rotation and shimmer were enabled. Optimization campaign and current-state update (supporting artifact removed).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base plus selective `HS_O3_FN` physics, shader, cull, lattice-distance/refinement and orientation kernels |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, source `7baf3cc430753e1a3693134150b282f212aa4e9e` (clean captured candidate) |
| Method | `HS_PROFILE`, 32-frame windows, 130 s capture, epoch stretched to 1200 revolutions (150 s). Exact runtime frames 2–444; setup frame 1 excluded. Scope/ISR summaries use complete windows 33–416. |
| Reproduce | `HS_TEENSY_PORT=COM4 HS_PROFILE_TREE=/c/work/Holosphere bash tools/profile_one.sh GSReactionDiffusion profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Single-effect image size: `FLASH: code:66352, data:252240, headers:9088   free for files:1703936` / `RAM1: variables:315136, code:28040, padding:4728   free for local variables:176384` / `RAM2: variables:520064  free for malloc/new:4224`.

Full-roster shipping Phantasm, separately built with profiling disabled: `FLASH: code:535552, data:746852, headers:8860   free for files:740352` / `RAM1: variables:314784, code:171928, padding:24680   free for local variables:12896` / `RAM2: variables:520064  free for malloc/new:4224`. The full-roster build passes the region-budget and layout gates; its RAM1 code is the shipping ITCM budget, separate from the single-effect comparison image.

Exactness cross-check: window frames 353–384, root 6,105,645,223 cycles ÷ 600 MHz versus measured wall sum 10,176,077 us, within **0.2 ppm**. Parser validation (supporting artifact removed) reports VALID.

## Frame cadence

README cells: peak 🔴 313.695, spilled 🔴 443/443 (100.00%).

**Runtime aggregate** (setup frame 1 excluded): mean render **262.576 ms/frame**, peak **313.695 ms** at frame 366, spilled **443/443 (100.00%)**. Runtime mean wall is 289.452 ms/frame, or 3.45 rendered frames/s across these live rows.

Startup setup render: **143.577 ms** at frame 1, before publication and excluded from runtime statistics.

A display window is 62.5 ms. The effect renders one 144×72 quadrant (10,368 positions), with four samples per shaded pixel. Every captured live frame exceeds the 16 fps render deadline. Peak render exceeds that deadline by 251.195 ms and requires 5.02× render speedup to fit one window. `canvas_buffer_wait` is synchronization idle until a display flip, not rendering work.

| Scope | Mean ms/frame | Worst window ms/frame |
|---|--:|--:|
| `grd_render` | 269.227 | 301.979 |
| `grd_simulate` | 106.832 | 107.017 |
| `grd_rasterize` | 159.123 | 191.963 |
| `grd_shader_draw` | 158.284 | 191.243 |
| `grd_cull_flags` | 0.466 | 0.781 |
| `grd_orient` | 0.373 | 0.374 |
| `canvas_buffer_wait` | 26.616 | 44.189 |

## Phase-by-phase readout

Phase schedule: the evolving reaction field changes active shading coverage. The capture has no lifecycle markers, so these windows describe measured workloads without assigning unobserved dissolve/reseed boundaries. The slower image advances fewer simulation frames within 130 seconds; this is a measured interval, not a lifetime worst-case bound.

### Early captured growth (frames 33–64)

```
frame                     226.29 ms 135.77 Mcyc 100.0%
  pov_preserve_half        135.8 us   81.5 kcyc   0.1% x1 135.8us/call
  grd_render              191.27 ms 114.76 Mcyc  84.5%
    grd_rasterize          80.99 ms  48.59 Mcyc  35.8%
      grd_shader_draw      79.84 ms  47.90 Mcyc  35.3% x1 79835.5us/call
      grd_cull_flags       780.8 us  468.5 kcyc   0.3% x1 780.8us/call
      grd_orient           373.6 us  224.1 kcyc   0.2% x1 373.6us/call
    grd_simulate          107.02 ms  64.21 Mcyc  47.3% x1 107017.0us/call
  rd_timeline_step          31.2 us   18.7 kcyc   0.0% x1 31.2us/call
  canvas_clear              85.0 us   51.0 kcyc   0.0% x1 85.0us/call
  canvas_buffer_wait       34.76 ms  20.86 Mcyc  15.4% x1 34764.4us/call
```

Wall min/avg/max = 155.566/226.287/253.163 ms. Shader work costs 79.835 ms/frame and simulation 107.017 ms/frame. Both contribute materially to the deadline overrun; higher active coverage raises the shader cost. These scope costs are window means, distinct from the exact individual-frame peak.

### Highest mean render window (frames 353–384)

```
frame                     318.00 ms 190.80 Mcyc 100.0%
  pov_preserve_half        126.9 us   76.2 kcyc   0.0% x1 126.9us/call
  grd_render              301.98 ms 181.19 Mcyc  95.0%
    grd_rasterize         191.96 ms 115.18 Mcyc  60.4%
      grd_shader_draw     191.24 ms 114.75 Mcyc  60.1% x1 191242.8us/call
      grd_cull_flags       346.1 us  207.7 kcyc   0.1% x1 346.1us/call
      grd_orient           373.5 us  224.1 kcyc   0.1% x1 373.5us/call
    grd_simulate          106.74 ms  64.05 Mcyc  33.6% x1 106744.3us/call
  rd_timeline_step          26.2 us   15.7 kcyc   0.0% x1 26.2us/call
  canvas_clear              85.0 us   51.0 kcyc   0.0% x1 85.0us/call
  canvas_buffer_wait       15.78 ms   9.47 Mcyc   5.0% x1 15784.7us/call
```

Wall min/avg/max = 293.012/318.002/374.591 ms. Shader work costs 191.243 ms/frame and simulation 106.744 ms/frame. Both contribute materially to the deadline overrun; higher active coverage raises the shader cost. These scope costs are window means, distinct from the exact individual-frame peak.

### Per-pixel figures

No `filter_blend` counter is present. The raster visits 10,368 quadrant positions, but active-pixel counts are unavailable. Dividing by the full quadrant would not measure cost per shaded pixel; no per-pixel cost is inferred.

## Column-ISR / DMA marshaling cost

Complete post-startup windows only; columns are calls/frame, per-call min/avg/max, and CPU share. Pack and submit are nested inside wake, so their shares are not added to its inclusive share.

```
isr_wake          5458.6/f 0.46/1.69/15.83 us  3.11% CPU
  isr_pack         682.3/f 6.23/7.07/11.04 us  1.63% CPU
  isr_dma_submit   682.3/f 0.62/0.94/8.28 us  0.22% CPU
```

- Pack averages 7.073 us/call; DMA submission averages 0.937 us/call. CPU shares use logged window elapsed time.
- SPI DMA transfers asynchronously. At 24 MHz, a 72-LED, 300-byte frame takes about 115 us including LPSPI framing; image plus trailing black frame takes about 230 us. These are wire-time calculations, not measured CPU submission times.
- Inclusive wake share 3.11% leaves about 60.557 ms of foreground CPU per 62.5 ms interval before other interrupts. Foreground scope times already include interrupts; the observed peak requires 5.02× render speedup to fit the wall deadline.

## Summary ranking

1. `grd_shader_draw` — 60.1% of the highest-render window, 191.243 ms/frame.
2. `grd_simulate` — 33.6% of the highest-render window, 106.744 ms/frame.
3. `grd_orient` — 0.1% of the highest-render window, 0.374 ms/frame.
4. `grd_cull_flags` — 0.1% of the highest-render window, 0.346 ms/frame.

Native quality/state evidence in the historical campaign predates the added shading features. It supplies no directly comparable native/WASM timing for this measured candidate.

## Caveats

- CYCCNT free-runs, so foreground scopes include ISR time; wake includes nested pack/submit work.
- `filter_blend` can inherit the first active parent and disappear with an inactive parent. This capture contains no such counter, and no per-pixel profiling scopes are enabled.
- Shipping selective-O3 covers GS physics, shading and hot-flag kernels plus shared lattice/orientation and renderer/driver paths. Global-O3 is a single-effect compiler comparison, not the full-roster shipping configuration.
- The extended epoch prevents reinitialization during capture. No dwell compression, simulation-speed override, or lifecycle completion is claimed.
- Both current images were captured on COM4 at source `7baf3cc430753e1a3693134150b282f212aa4e9e` with an empty source diff. The candidate includes the noise modifier speed and pigment scratch-lifetime corrections. The prior early-morning captures used COM3 and earlier artwork; their timing delta is not a controlled optimization comparison.
- Setup frame 1 is excluded from runtime rows; the entire first window is excluded from scope/ISR summaries. Complete individual frame rows in the unfinished final window remain in runtime statistics. Integer-microsecond ISR totals introduce small quantization error.
- Portable evidence manifest (supporting artifact removed) retains capture, compiler/build/environment records and original hashes. Host process environment dictionaries are removed; footer environment hashes identify the sanitized retained dumps, while measurement rows are unchanged. ELF/map artifacts remain in the provenance-named local archive.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1200`. `just profile GSReactionDiffusion` is the basic shortcut; use the Setup command for this duration and epoch under the shared-device lock.

Source reachability: capture `7baf3cc4307` maps to landed `daf63d812` on a
different base (319 files differ); the captured tree is not available from the
published branch. Supporting artifacts are no longer retained.
