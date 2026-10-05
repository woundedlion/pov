# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-10-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).

Raw capture: `C:/work/temp/gs-perf-20261004/final-shipping.log`, captured 2026-10-05 00:38 local time on COM4. Replaces the 2026-09-30 GSReactionDiffusion report. Build, compiler, environment and ELF hashes are retained beside the capture and in its provenance footer.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env, `-Os` base with selective `HS_O3_FN` chemistry, pigment, palette/shader, cull, lattice/refinement and orientation paths; no per-pixel profiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, clean source `c25f15b0827590db54f8f636fbc0c7efccbeb89b` |
| Method | `HS_PROFILE`, 32-frame windows, 130 s capture; epoch = 1200 revolutions (150 s). Exact runtime frames 2–2033; frame 1 excluded. Scope/ISR summaries use complete windows 33–2016. |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Single-effect profile image: `FLASH: code:70388, data:344916, headers:8628   free for files:1607684` / `RAM1: variables:315328, code:36792, padding:28744   free for local variables:143424` / `RAM2: variables:520064  free for malloc/new:4224`.

Full-roster shipping Phantasm, separately built with profiling disabled: `FLASH: code:565620, data:845828, headers:8836   free for files:611332` / `RAM1: variables:314784, code:180008, padding:16600   free for local variables:12896` / `RAM2: variables:520064  free for malloc/new:4224`.

Compared with the current-session baseline Phantasm, flash code grows 10,192 bytes and flash data 92,316 bytes. RAM1 code grows 8,816 bytes, offset by less padding; RAM1 variables and free stack space remain 314,784 and 12,896 bytes. RAM2 remains 520,064 bytes used / 4,224 bytes free. The profile-only image has 143,424 bytes free RAM1; this is separate from the full-roster shipping budget.

| Image | Flash code before → after | Flash data before → after | RAM1 code before → after | RAM1 free before → after |
|---|--:|--:|--:|--:|
| Single-effect profile | 66,556 → 70,388 | 253,220 → 344,916 | 28,168 → 36,792 | 176,384 → 143,424 |
| Full-roster Phantasm | 555,428 → 565,620 | 753,512 → 845,828 | 171,192 → 180,008 | 12,896 → 12,896 |

Exactness cross-check: frames 385–416, root 1,198,505,901 cycles ÷ 600 MHz versus wall sum 1,997,516 us, within **3.1 ppm**. `final-validation.txt` reports VALID, including complete per-frame telemetry and no epoch reset.

## Frame cadence

README cells: peak 🟢 39.271, spilled 🟢 0/2032 (0.00%).

**Runtime aggregate:** mean render **32.758 ms/frame**, peak **39.271 ms** at frame **395**, spilled **0/2032 (0%)**. All 2032 live frames are below 40 ms; measured peak margin is **0.729 ms**. Mean wall is 62.440 ms/frame. Frame 1 render is 21.593 ms and is excluded from runtime statistics; initialization scopes logged before frame publication are excluded as well.

A display window is 62.5 ms at 480 RPM. Each frame renders one 144×72 quadrant (10,368 positions), retaining four geometric coverage samples per shaded pixel. The measured peak leaves 23.229 ms before the display deadline. `canvas_buffer_wait` is synchronization idle, not rendering work. Individual wall times vary around the cadence as rendering costs and display phase change.

The current-session baseline peak was **311.041 ms**, with 443/443 live frames spilling; this candidate reduces measured peak by **87.37%** (7.92× faster). The superseded September 30 report recorded 313.695 ms. These are captured peaks, not bounds over every parameter setting or random seed.

Complete-window scope averages cover 1984 frames, excluding the first 32-frame window; they differ from the live-row aggregate above.

| Scope | Mean ms/frame | Worst window ms/frame |
|---|--:|--:|
| `grd_render` | 32.666 | 38.505 |
| `grd_simulate` | 9.973 | 10.225 |
| `grd_physics` | 5.107 | 5.123 |
| `grd_pigment` | 4.118 | 4.362 |
| `grd_rasterize` | 19.808 | 25.059 |
| `grd_shader_draw` | 17.424 | 22.864 |
| `grd_cull_flags` | 0.525 | 0.937 |
| `grd_orient` | 1.859 | 1.868 |
| `grd_color_noise` | 2.711 | 2.733 |
| `canvas_buffer_wait` | 29.537 | 40.206 |

## Phase-by-phase readout

The capture includes growth, denser shading and four windows with runtime seed/palette scopes. It has no explicit lifecycle markers; the window descriptions below follow measured workload and scope activity, without assigning an exact dissolve boundary from timing alone. Trees preserve logged nesting and express each counter as a per-frame mean; percentages use the frame root.

### Early captured growth (frames 33–64)

```
frame                        62.596 ms 37.557 Mcyc 100.0%
  pov_preserve_half           0.145 ms  0.087 Mcyc   0.2%
  grd_render                 25.381 ms 15.229 Mcyc  40.5%
    grd_rasterize            13.156 ms  7.893 Mcyc  21.0%
      grd_shader_draw        10.503 ms  6.302 Mcyc  16.8%
      grd_cull_flags          0.794 ms  0.476 Mcyc   1.3%
      grd_orient              1.858 ms  1.115 Mcyc   3.0%
    grd_simulate              9.528 ms  5.717 Mcyc  15.2%
      grd_physics             5.121 ms  3.072 Mcyc   8.2%
      grd_pigment             3.661 ms  2.196 Mcyc   5.8%
    grd_color_noise           2.697 ms  1.618 Mcyc   4.3%
  rd_timeline_step            0.033 ms  0.020 Mcyc   0.1%
  canvas_clear                0.084 ms  0.050 Mcyc   0.1%
  canvas_buffer_wait         36.951 ms 22.170 Mcyc  59.0%
```

Wall min/avg/max = 58.303/62.596/66.955 ms. The shader is cheaper at this early workload. Chemistry still runs six substeps per simulation frame; pigment transport runs once per simulation frame.

### Highest mean render window (frames 385–416)

```
frame                        62.422 ms 37.453 Mcyc 100.0%
  pov_preserve_half           0.140 ms  0.084 Mcyc   0.2%
  grd_render                 38.505 ms 23.103 Mcyc  61.7%
    grd_rasterize            25.059 ms 15.035 Mcyc  40.1%
      grd_shader_draw        22.864 ms 13.719 Mcyc  36.6%
      grd_cull_flags          0.335 ms  0.201 Mcyc   0.5%
      grd_orient              1.859 ms  1.116 Mcyc   3.0%
    grd_simulate             10.188 ms  6.113 Mcyc  16.3%
      grd_physics             5.117 ms  3.070 Mcyc   8.2%
      grd_pigment             4.319 ms  2.592 Mcyc   6.9%
    grd_color_noise           2.722 ms  1.633 Mcyc   4.4%
  rd_timeline_step            0.026 ms  0.015 Mcyc   0.0%
  canvas_clear                0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait         23.666 ms 14.199 Mcyc  37.9%
```

Wall min/avg/max = 61.467/62.422/63.451 ms. This window also contains the exact 39.271 ms peak at frame 395. Shader work is the largest cost; the idle buffer wait keeps the display synchronized. Physics records six calls/frame at 852.795 us/call; pigment records one call/frame at 4.319 ms/call.

### First runtime reseed window (frames 449–480)

```
frame                        62.290 ms 37.374 Mcyc 100.0%
  pov_preserve_half           0.143 ms  0.086 Mcyc   0.2%
  grd_render                 23.328 ms 13.997 Mcyc  37.5%
    grd_rasterize            10.636 ms  6.382 Mcyc  17.1%
      grd_shader_draw         7.920 ms  4.752 Mcyc  12.7%
      grd_cull_flags          0.914 ms  0.548 Mcyc   1.5%
      grd_orient              1.802 ms  1.081 Mcyc   2.9%
    grd_simulate              8.800 ms  5.280 Mcyc  14.1%
      grd_physics             4.960 ms  2.976 Mcyc   8.0%
      grd_pigment             3.116 ms  1.869 Mcyc   5.0%
    grd_color_noise           2.710 ms  1.626 Mcyc   4.4%
  rd_timeline_step            0.032 ms  0.019 Mcyc   0.1%
  canvas_clear                0.084 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait         38.701 ms 23.221 Mcyc  62.1%
grd_color_palette             0.412 ms  0.247 Mcyc   0.7%
grd_seed_reaction             0.557 ms  0.334 Mcyc   0.9%
```

Wall min/avg/max = 60.208/62.290/66.321 ms. The 31 raster/simulation calls across 32 frames reflect the staged seed lifecycle. `grd_seed_reaction` and `grd_color_palette` are tagged MIXED-PARENT: their top-level printed totals are attribution diagnostics, not additional exclusive work to add to `grd_render`.

All runtime reseed windows remain below 40 ms:

| Frames | Seed-scope calls | Palette-scope calls | `grd_render` mean ms | Peak frame render ms |
|---|--:|--:|--:|--:|
| 449–480 | 2 | 6 | 23.328 | 29.918 |
| 865–896 | 2 | 6 | 21.831 | 27.021 |
| 1281–1312 | 2 | 6 | 25.420 | 33.696 |
| 1793–1824 | 2 | 6 | 22.787 | 30.512 |

### Per-pixel figures

The raster visits 10,368 quadrant positions. There is no active-pixel or `filter_blend` counter in this capture, so no cost per shaded pixel is inferred. Per-pixel diagnostic scopes are disabled.

## Column-ISR / DMA marshaling cost

Complete post-startup windows only. Values are calls/frame, per-call min/weighted mean/max, and CPU share over logged window elapsed time. Pack and submit nest inside wake; their shares are not added to the inclusive wake share.

```
isr_wake        1152.2/f  0.54/1.66/19.17 us  3.07% CPU
  isr_pack         144.0/f  6.24/7.10/12.29 us  1.64% CPU
  isr_dma_submit   144.0/f  0.58/0.94/11.79 us  0.22% CPU
```

Pack averages 7.100 us/call; DMA submission averages 0.940 us/call. Submission measures CPU setup, while SPI wire transfer proceeds asynchronously and is not timed by this scope. Inclusive wake share 3.07% corresponds to about 60.583 ms of foreground CPU per 62.5 ms interval before other interrupts. Foreground scope measurements already include those ISR interruptions.

## Summary ranking

1. `grd_shader_draw` — 36.6% of the highest-render window, 22.864 ms/frame.
2. `grd_physics` — 8.2% of the highest-render window, 5.117 ms/frame.
3. `grd_pigment` — 6.9% of the highest-render window, 4.319 ms/frame.
4. `grd_color_noise` — 4.4% of the highest-render window, 2.722 ms/frame.
5. `grd_orient` — 3.0% of the highest-render window, 1.859 ms/frame.

The regression's main costs were repeated nonlinear palette transforms inside the four-sample shader and pigment transport repeated for every chemistry substep. The optimized path caches palette colors, evaluates pixel color once, and advances pigment once per frame. Chemistry and raster geometry retain their resolution. No directly comparable WASM/native timing is claimed.

## Caveats

- This establishes the under-40 ms result for the default controls and captured lifecycle on COM4, with 0.729 ms measured margin. It does not establish a universal peak over all controls, seeds, temperatures or devices. Larger hue/shimmer ranges use the procedural color fallback and may cost more.
- Pigment flow intentionally changes: one frozen-neighbor transport update per frame replaces six transport updates. All six chemistry substeps remain. Shading uses nearest-node pigment, a 16×15 float RGB lookup per palette, one color evaluation from mean concentration, and paired reciprocal estimates for concentration. All four geometric coverage samples remain. These are measured visual approximations, not bit-identical color output to the original effect.
- The seed transition publishes one black frame while staging the next reaction. The following frame completes seeding; palette rows rebuild incrementally with procedural fallback for rows not ready. Initialization is outside the runtime limit.
- CYCCNT free-runs, so foreground scopes include ISR time. `grd_seed_reaction` and `grd_color_palette` are MIXED-PARENT and must not be summed as exclusive siblings.
- The build keeps shipping selective-O3 behavior; no global-O3 comparison or per-pixel instrumentation is included. The extended epoch avoids harness teardown during capture; no simulation-speed override or dwell compression was used.
- Frame 1 and the first counter window are excluded from the corresponding runtime/scope summaries. Seventeen complete raw frame rows after the last counter window remain in runtime statistics. ISR total microseconds have integer quantization.
- Source was clean at `c25f15b08`. Raw log, provenance, build logs, environment records and ELF archive remain in `C:/work/temp/gs-perf-20261004`; `final-shipping-metrics.json` records the capture SHA-256 and computed statistics.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=GSReactionDiffusion`, `HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1200`. `just profile GSReactionDiffusion` is the basic shortcut; use the Setup command for this duration and epoch under the shared-device lock.
