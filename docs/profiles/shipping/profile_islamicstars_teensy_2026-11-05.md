# IslamicStars on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile IslamicStars`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/islamicstars_ship.log`. Captured 2026-11-04 23:39 America/Los_Angeles on COM4. Replaces profile_islamicstars_teensy_2026-10-30.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | IslamicStars, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 306 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh IslamicStars profile 306 16` |

Image size:

```text
teensy_size:   FLASH: code:190792, data:197628, headers:8892   free for files:1634304
teensy_size:    RAM1: variables:314688, code:30856, padding:1912   free for local variables:176832
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `936c7b6bf6181998aceed3cdc9f79c8cdecbed9c58e79376a0af1888e3826103`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 2049–2064, root cyc/600 versus measured wall sum differs by **15.182 ppm**. This fails the 5 ppm check; see validator diagnostics below.

Hot-layout gate in the measured profile build:

```text
  ok   IslamicStars             worst 2 lines/set (0 sets over)  gated 20,474 B  overflow 45,900 B  full-set worst 5
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **59.311 ms** (frame 2056), **0/4863 spilled** (0.00%). Setup frame 1: 16.664 ms. Complete windows cover frames 1–4864; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 24.753 ms/frame. `is_timeline_step` averages 24.456 ms/frame over complete windows; its worst window is 56.511 ms/frame (frames 2049–2064).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 8.624–59.311 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     56.323 ms 33.794 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x0.94/frame; 147.586 us/call
  is_timeline_step                        10.247 ms  6.148 Mcyc  18.2%
    is_draw_shape                         10.069 ms  6.041 Mcyc  17.9%
      is_mesh_scan                         9.833 ms  5.900 Mcyc  17.5%
        scan_mesh_raster                   9.370 ms  5.622 Mcyc  16.6%
          x7.94/frame; 1180.524 us/call
        scan_face_setup                    0.221 ms  0.133 Mcyc   0.4%
          x7.94/frame; 27.890 us/call
      is_face_offsets                      0.013 ms  0.008 Mcyc   0.0%
        x1.00/frame; 13.179 us/call
      is_mesh_transform                    0.018 ms  0.011 Mcyc   0.0%
        x1.00/frame; 17.675 us/call
      mesh_visibility                      0.020 ms  0.012 Mcyc   0.0%
        x1.00/frame; 19.584 us/call
  is_ripple_prepare                        0.003 ms  0.002 Mcyc   0.0%
    x1.00/frame; 2.804 us/call
  canvas_clear                             0.104 ms  0.062 Mcyc   0.2%
    x1.00/frame; 103.613 us/call
  canvas_buffer_wait                      45.808 ms 27.485 Mcyc  81.3%
    x1.00/frame; 45807.631 us/call
```

Wall min/avg/max: 9.389/56.322/63.658 ms. Window mean render is 10.515 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 2049–2064)

```text
scan_face_setup                            0.400 ms  0.240 Mcyc   0.6%
  x17.88/frame; 22.356 us/call
mesh_copy_capture                          0.345 ms  0.207 Mcyc   0.6%
  x22.19/frame; 15.559 us/call
mesh_capture                               2.814 ms  1.688 Mcyc   4.5%
  x5.00/frame; 562.718 us/call
mesh_classify                              1.484 ms  0.890 Mcyc   2.4%
  x1.00/frame; 1484.123 us/call
frame                                     62.226 ms 37.335 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.467 us/call
  is_timeline_step                        56.511 ms 33.907 Mcyc  90.8%
    is_draw_shape                         56.300 ms 33.780 Mcyc  90.5%
      is_mesh_scan                        46.566 ms 27.940 Mcyc  74.8%
        scan_face_sweep                    5.160 ms  3.096 Mcyc   8.3%
          scan_sweep_column                1.953 ms  1.172 Mcyc   3.1%
            x143.25/frame; 13.636 us/call
        mesh_sweep_paths                   0.126 ms  0.075 Mcyc   0.2%
          x1.00/frame; 125.648 us/call
        scan_mesh_raster                  26.773 ms 16.064 Mcyc  43.0%
          x402.19/frame; 66.568 us/call
        scan_face_setup                    5.785 ms  3.471 Mcyc   9.3%
          x426.88/frame; 13.551 us/call
      is_face_offsets                      0.035 ms  0.021 Mcyc   0.1%
        x1.00/frame; 35.257 us/call
      is_mesh_transform                    2.576 ms  1.546 Mcyc   4.1%
        x1.00/frame; 2576.158 us/call
      mesh_visibility                      1.226 ms  0.736 Mcyc   2.0%
        x1.00/frame; 1226.300 us/call
  is_ripple_prepare                        0.020 ms  0.012 Mcyc   0.0%
    x1.00/frame; 19.748 us/call
  canvas_clear                             0.094 ms  0.056 Mcyc   0.2%
    x1.00/frame; 93.559 us/call
  canvas_buffer_wait                       5.407 ms  3.244 Mcyc   8.7%
    x1.00/frame; 5407.316 us/call
```

Wall min/avg/max: 59.414/62.224/65.371 ms. Window mean render is 56.818 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

IslamicStars clean holds require `scan_mesh_raster` calls = window frames × finished F from the Built Shape log. Spawning Shape counts tagged seed describe the initial mesh and are labeled separately. Rows lacking that full-face-count observation are marked as anchor estimates, not clean holds.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 seed V=20 E=30 F=12 I=60; built V=3240 E=4320 F=1082 I=8640 | 22/0 | unavailable | 56.511 | 56.818 | 16.07 | anchor; no full-face-count hold |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 seed V=60 E=90 F=32 I=180; built V=1620 E=2160 F=542 I=4320 | 18/0 | unavailable | 46.751 | 47.056 | 15.98 | anchor; no full-face-count hold |
| truncatedIcosidodecahedron_truncate50d_ambo_dual seed V=120 E=180 F=62 I=360; built V=542 E=1080 F=540 I=2160 | 19/0 | unavailable | 44.321 | 45.288 | 15.92 | anchor; no full-face-count hold |
| truncatedCube_hk90_meta seed V=24 E=36 F=14 I=72; built V=218 E=648 F=432 I=1296 | 23/0 | unavailable | 42.459 | 42.747 | 15.90 | anchor; no full-face-count hold |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 seed V=12 E=30 F=20 I=60; built V=1620 E=2160 F=542 I=4320 | 18/0 | unavailable | 39.665 | 39.962 | 15.99 | anchor; no full-face-count hold |
| icosahedron_snub_relax_truncate033_hankin62 seed V=12 E=30 F=20 I=60; built V=1350 E=1800 F=452 I=3600 | 18/0 | unavailable | 32.759 | 33.055 | 16.01 | anchor; no full-face-count hold |
| snubDodecahedron_truncate5d_ambo_dual seed V=60 E=150 F=92 I=300; built V=452 E=900 F=450 I=1800 | 19/0 | unavailable | 32.276 | 32.556 | 16.02 | anchor; no full-face-count hold |
| truncatedIcosahedron_hk58_chamfer63 seed V=60 E=90 F=32 I=180; built V=990 E=1440 F=452 I=2880 | 15/0 | unavailable | 31.931 | 32.224 | 16.05 | anchor; no full-face-count hold |
| dodecahedron_hk72_ambo_dual_hk20 seed V=20 E=30 F=12 I=60; built V=540 E=720 F=182 I=1440 | 22/0 | unavailable | 31.115 | 31.421 | 16.03 | anchor; no full-face-count hold |
| dodecahedron_bevel2_relax_gyro seed V=20 E=30 F=12 I=60; built V=542 E=900 F=360 I=1800 | 22/0 | unavailable | 30.714 | 30.994 | 16.00 | anchor; no full-face-count hold |
| rhombicuboctahedron_hk63_ambo_hk63 seed V=24 E=48 F=26 I=96; built V=864 E=1152 F=290 I=2304 | 18/0 | unavailable | 30.157 | 30.453 | 16.04 | anchor; no full-face-count hold |
| icosahedron_ambo_truncate033_hankin59 seed V=12 E=30 F=20 I=60; built V=540 E=720 F=182 I=1440 | 17/0 | unavailable | 27.814 | 28.126 | 16.05 | anchor; no full-face-count hold |
| dodecahedron_hk54_ambo_hk72 seed V=20 E=30 F=12 I=60; built V=540 E=720 F=182 I=1440 | 17/0 | unavailable | 26.101 | 26.399 | 16.06 | anchor; no full-face-count hold |
| icosidodecahedron_truncate5d_ambo_dual seed V=30 E=60 F=32 I=120; built V=182 E=360 F=180 I=720 | 20/0 | unavailable | 25.585 | 25.867 | 15.99 | anchor; no full-face-count hold |
| dodecahedron_hk62_ambo_hk62 seed V=20 E=30 F=12 I=60; built V=540 E=720 F=182 I=1440 | 18/0 | unavailable | 25.412 | 25.707 | 16.03 | anchor; no full-face-count hold |
| octahedron_hk17_ambo_hk73 seed V=6 E=12 F=8 I=24; built V=216 E=288 F=74 I=576 | 18/0 | unavailable | 24.927 | 25.237 | 16.02 | anchor; no full-face-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| dodecahedron_hk35_ambo_hk62_ambo_relax_hk42 | 59.311 | 0/352 | 16.00 |
| truncatedIcosidodecahedron_truncate50d_ambo_dual | 50.852 | 0/312 | 16.00 |
| truncatedCube_hk90_meta | 50.306 | 0/368 | 16.00 |
| truncatedIcosahedron_ambo_relax_truncate001_hankin59 | 49.682 | 0/288 | 16.00 |
| icosahedron_truncate25_hk60_ambo_relax100_hk60 | 41.826 | 0/289 | 16.00 |
| octahedron_hk17_ambo_hk73 | 35.770 | 0/280 | 16.00 |
| snubDodecahedron_truncate5d_ambo_dual | 34.811 | 0/312 | 16.00 |
| icosahedron_snub_relax_truncate033_hankin62 | 34.753 | 0/288 | 16.00 |
| dodecahedron_hk72_ambo_dual_hk20 | 34.635 | 0/352 | 16.00 |
| dodecahedron_bevel2_relax_gyro | 34.191 | 0/352 | 16.00 |
| truncatedIcosahedron_hk58_chamfer63 | 33.703 | 0/248 | 16.00 |
| icosahedron_ambo_truncate033_hankin59 | 32.206 | 0/272 | 16.00 |
| rhombicuboctahedron_hk63_ambo_hk63 | 32.158 | 0/280 | 16.00 |
| dodecahedron_hk54_ambo_hk72 | 27.454 | 0/280 | 16.00 |
| icosidodecahedron_truncate5d_ambo_dual | 27.083 | 0/312 | 16.00 |
| dodecahedron_hk62_ambo_hk62 | 26.361 | 0/278 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.12/frame
  min/avg/max 0.516/1.980/15.585 us; cpu 3.65%
  isr_pack         143.88/frame
    min/avg/max 8.773/9.331/11.986 us; cpu 2.14%
  isr_dma_submit   143.88/frame
    min/avg/max 0.616/0.936/2.871 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=999217 us
isr_wake       n=18418 cyc min/avg/max=310/1188/9351 ns min/avg/max=516/1980/15585 total=36494 us cpu=3.65%
isr_pack       n=2302 cyc min/avg/max=5264/5599/7192 ns min/avg/max=8773/9331/11986 total=21482 us cpu=2.14%
isr_dma_submit n=2302 cyc min/avg/max=370/562/1723 ns min/avg/max=616/936/2871 total=2156 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.65% leaves approximately 60.219 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `is_timeline_step` — 56.511 ms/frame, 90.82% of this root window.
2. `is_draw_shape` — 56.300 ms/frame, 90.48% of this root window.
3. `is_mesh_scan` — 46.566 ms/frame, 74.83% of this root window.
4. `scan_mesh_raster` — 26.773 ms/frame, 43.03% of this root window.
5. `scan_face_setup` — 5.785 ms/frame, 9.30% of this root window.

README cells: peak 🟢 59.311 (16), spilled 🟢 0/4863 (0.00%).

README setup frame 1: 16.664 ms. Captured 2026-11-04 23:39 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under its first caller and can be hidden when that parent has no calls; per-pixel instrumentation adds overhead when enabled.
- Selective `HS_O3` regions are active through this effect’s instrumented shipping draw path; exact flags and compiler/ABI attestation are adjacent to the raw capture.
- No dwell compression, transition-speed override, ordered-cycle override or fixed preset is added to the stock sweep. The default one-hour epoch is used.
- Captured from a clean detached source worktree. Generated build artifacts are untracked outputs.
- IslamicStars uses default transition speed, matching 11-04; the 10-30 TRANS_SPEED=4 report has a different schedule.

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
  [FAIL] root cyc/600 vs wall sum within 5 ppm (15.2 ppm, frames 2049-2064)

INVALID: IslamicStars
```

Validation flags: [FAIL] shape markers advance and return to the first (cycle closed): 16 distinct; [FAIL] root cyc/600 vs wall sum within 5 ppm (15.2 ppm, frames 2049-2064). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=IslamicStars`, `HS_PROFILE_WINDOW=16`. `just profile IslamicStars` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
