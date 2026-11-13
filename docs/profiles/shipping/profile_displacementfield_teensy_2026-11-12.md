# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-11-12, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `C:/work/temp/df-core-device/final/fullcycle.log`. Captured 2026-11-12 11:23 America/Los_Angeles on COM4 (serial 6362390). Replaces `profile_displacementfield_teensy_2026-11-07.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, newlib-nano and selective `HS_O3`; ball-field evaluators and distorted-ring stack kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, clean source `182eb03a5c86145197a3723386ab0f0c66f12a3e` |
| Method | One 190 s full-cycle pass, 16-frame counter windows, normal preset timing; `HS_PROFILE_EPOCH_REVS=2400` (300 s) |
| Reproduce | `HS_TEENSY_PORT=COM4 HS_PROFILE_TREE=C:/work/Holosphere bash tools/profile_one.sh DisplacementField profile 190 16 -D HS_PROFILE_EPOCH_REVS=2400` |

Image size:

```text
teensy_size:   FLASH: code:116760, data:150460, headers:8236   free for files:1756160
teensy_size:    RAM1: variables:314464, code:25384, padding:7384   free for local variables:177056
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `315d08d58f7483105ff028722e772abf63bad8375b18b3061abcc7df657457f8`. Phantasm ELF SHA-256: `791cb279df6d14c6754f8e551a5135000f28398973754b7ca768ae368dd015f0`. Artifacts: `C:/work/temp/df-core-device/final/artifacts/fullcycle_182eb03a5c86_315d08d58f74`; adjacent provenance, envdump, build logs and source status retain compiler/ABI and build inputs.

Exactness cross-check: frames 2113–2128, root 598,269,370 cycles ÷ 600 MHz = 997115.617 µs versus measured wall sum 997,116 µs: **0.384 ppm**, below the 5 ppm gate.

## Frame cadence

**Live render, setup frame 1 excluded:** peak **45.164 ms** (frame 603), **0/3015 spilled** (0.00%). Setup frame 1: **18.393 ms**. Every complete per-frame record is included, including the trailing partial counter window; live transitions remain included.

Complete counter windows cover 3008 frames across 188 windows. Mean live render is 36.674 ms/frame. The worst `df_timeline_step` window averages 44.164 ms/frame (frames 2113–2128).

At 480 RPM each half revolution is 62.5 ms, yielding 16 fps when rendering fits one window. The effect renders one quadrant, 144×72 = 10,368 pixels. Noise, balls and transitions all fit this cadence. `canvas_buffer_wait` is the round-up idle to the next display flip.

## Phase-by-phase readout

Normal scheduling fades noise in, holds and fades it out, spawns and drains balls, then returns to noise. The capture confirms noise → balls → noise. Counter trees preserve raw nesting and duplicate counter instances; mixed-parent/duplicate-name rows are inclusive and must not be summed as exclusive work.

### Startup / noise entry (frames 1–16)

```text
frame                        57.775 ms  34.665 Mcyc 100.0%
  pov_preserve_half         138.605 us  83.163 kcyc   0.2%
    x0.94/frame; 147.846 us/call
  df_timeline_step           25.184 ms  15.111 Mcyc  43.6%
    df_draw_rings            25.076 ms  15.046 Mcyc  43.4%
      df_fused_scan           9.865 ms   5.919 Mcyc  17.1%
        x0.56/frame; 17536.902 us/call
      df_lut_bake             7.182 ms   4.309 Mcyc  12.4%
        df_hue_table_prep   218.224 us 130.934 kcyc   0.4%
          x7.94/frame; 27.493 us/call
        df_octave_noise       4.931 ms   2.959 Mcyc   8.5%
          x26.81/frame; 183.922 us/call
      df_chunk_cull         445.532 us 267.319 kcyc   0.8%
        x26.88/frame; 16.578 us/call
      df_lut_bake             1.144 ms 686.255 kcyc   2.0%
        df_hue_table_prep    10.829 us   6.497 kcyc   0.0%
          x5.25/frame; 2.063 us/call
      df_chunk_cull         273.179 us 163.907 kcyc   0.5%
        x18.00/frame; 15.177 us/call
      df_fused_scan           4.901 ms   2.940 Mcyc   8.5%
        ring_stack_table      1.604 ms 962.419 kcyc   2.8%
          x1.50/frame; 1069.355 us/call
        ring_stack_table     87.132 us  52.279 kcyc   0.2%
          x0.06/frame; 1394.115 us/call
    animation_noise_clock     0.169 us    101.5 cyc   0.0%
      x1.00/frame; 0.169 us/call
  canvas_clear               99.571 us  59.743 kcyc   0.2%
    x1.00/frame; 99.571 us/call
  canvas_buffer_wait         32.324 ms  19.395 Mcyc  55.9%
    x1.00/frame; 32324.418 us/call
  df_prepare_fields           0.722 us    433.1 cyc   0.0%
    x1.00/frame; 0.722 us/call
```

Wall min/avg/max = 15.512/57.774/79.958 ms. Mean render is 25.450 ms/frame. Display synchronization accounts for the wait; wall maxima are not render spills.

### Held noise, worst timeline window (frames 2113–2128)

```text
ring_stack_table              2.123 ms   1.274 Mcyc   3.4%
  x2.00/frame; 1061.544 us/call
frame                        62.320 ms  37.392 Mcyc 100.0%
  pov_preserve_half         142.548 us  85.529 kcyc   0.2%
    x1.00/frame; 142.548 us/call
  df_timeline_step           44.164 ms  26.498 Mcyc  70.9%
    df_draw_rings            44.027 ms  26.416 Mcyc  70.6%
      df_fused_scan          26.031 ms  15.619 Mcyc  41.8%
        x1.00/frame; 26031.163 us/call
      df_lut_bake            15.652 ms   9.391 Mcyc  25.1%
        df_hue_table_prep     1.219 ms 731.169 kcyc   2.0%
          x23.56/frame; 51.718 us/call
        df_octave_noise      10.564 ms   6.339 Mcyc  17.0%
          x48.00/frame; 220.090 us/call
      df_chunk_cull         833.722 us 500.233 kcyc   1.3%
        x48.00/frame; 17.369 us/call
    animation_noise_clock     0.169 us    101.5 cyc   0.0%
      x1.00/frame; 0.169 us/call
  canvas_clear               93.486 us  56.092 kcyc   0.2%
    x1.00/frame; 93.486 us/call
  canvas_buffer_wait         17.887 ms  10.732 Mcyc  28.7%
    x1.00/frame; 17887.390 us/call
  df_prepare_fields           0.772 us    463.5 cyc   0.0%
    x1.00/frame; 0.772 us/call
```

Wall min/avg/max = 61.410/62.319/63.508 ms. Mean render is 44.432 ms/frame. Display synchronization accounts for the wait; wall maxima are not render spills.

### Balls, worst timeline window (frames 2977–2992)

```text
ring_stack_table              1.852 ms   1.111 Mcyc   3.0%
  x2.00/frame; 926.025 us/call
frame                        62.246 ms  37.348 Mcyc 100.0%
  pov_preserve_half         143.474 us  86.085 kcyc   0.2%
    x1.00/frame; 143.474 us/call
  df_timeline_step           36.125 ms  21.675 Mcyc  58.0%
    df_draw_rings            35.892 ms  21.535 Mcyc  57.7%
      df_fused_scan          20.667 ms  12.400 Mcyc  33.2%
        x1.00/frame; 20666.756 us/call
      df_lut_bake            12.996 ms   7.798 Mcyc  20.9%
        df_hue_table_prep   334.597 us 200.758 kcyc   0.5%
          x37.44/frame; 8.937 us/call
      df_chunk_cull         638.909 us 383.345 kcyc   1.0%
        x44.38/frame; 14.398 us/call
    animation_noise_clock     0.229 us    137.2 cyc   0.0%
      x1.00/frame; 0.229 us/call
  canvas_clear               92.321 us  55.393 kcyc   0.1%
    x1.00/frame; 92.321 us/call
  canvas_buffer_wait         25.827 ms  15.496 Mcyc  41.5%
    x1.00/frame; 25826.670 us/call
  df_prepare_fields          19.137 us  11.482 kcyc   0.0%
    x1.00/frame; 19.137 us/call
```

Wall min/avg/max = 60.606/62.246/63.659 ms. Mean render is 36.420 ms/frame. Display synchronization accounts for the wait; wall maxima are not render spills.

### Preset transition (frames 2801–2816)

```text
frame                        61.172 ms  36.703 Mcyc 100.0%
  pov_preserve_half         146.454 us  87.872 kcyc   0.2%
    x1.00/frame; 146.454 us/call
  df_timeline_step           18.145 ms  10.887 Mcyc  29.7%
    df_draw_rings            18.014 ms  10.808 Mcyc  29.4%
      df_fused_scan           8.392 ms   5.035 Mcyc  13.7%
        x0.62/frame; 13427.106 us/call
      df_lut_bake             3.677 ms   2.206 Mcyc   6.0%
        df_hue_table_prep   178.004 us 106.802 kcyc   0.3%
          x5.81/frame; 30.624 us/call
        df_octave_noise       2.509 ms   1.505 Mcyc   4.1%
          x11.06/frame; 226.808 us/call
      df_chunk_cull         228.479 us 137.088 kcyc   0.4%
        x13.06/frame; 17.491 us/call
      df_lut_bake           828.883 us 497.330 kcyc   1.4%
        df_hue_table_prep    10.749 us   6.449 kcyc   0.0%
          x5.38/frame; 2.000 us/call
      df_chunk_cull         194.649 us 116.789 kcyc   0.3%
        x12.25/frame; 15.890 us/call
      df_fused_scan           3.677 ms   2.206 Mcyc   6.0%
        ring_stack_table      1.190 ms 714.253 kcyc   1.9%
          x1.62/frame; 732.567 us/call
    animation_noise_clock     0.170 us    102.2 cyc   0.0%
      x1.00/frame; 0.170 us/call
  canvas_clear               90.805 us  54.483 kcyc   0.1%
    x1.00/frame; 90.805 us/call
  canvas_buffer_wait         42.747 ms  25.648 Mcyc  69.9%
    x1.00/frame; 42746.920 us/call
  df_prepare_fields           1.637 us    982.5 cyc   0.0%
    x1.00/frame; 1.637 us/call
```

Wall min/avg/max = 46.815/61.171/63.027 ms. Mean render is 18.425 ms/frame. Display synchronization accounts for the wait; wall maxima are not render spills.

### Per-preset table

Worst clean complete timeline window per preset, excluding setup and windows that straddle an advance. Blend counting is disabled in the shipping timing image. Both presets were visited and the capture returned to noise.

| # | Preset | Clean windows | Blended px/frame | Timeline ms | Render ms | fps |
|---|---|---:|---|---:|---:|---:|
| 1 | Noise | 55 | not captured | 44.164 | 44.432 | 16 |
| 2 | Balls | 74 | not captured | 36.125 | 36.420 | 16 |

### Per-pixel figures

Per-pixel blend scopes are disabled. Blended pixels/frame, coverage, cycles/blend and scan cycles/blended pixel are unavailable; quadrant area is not substituted for measured blended work.

## Column-ISR / DMA marshaling cost

Last complete window, 16 frames; separately measured ISR interval 997,995 µs:

```text
isr_wake        1149.75/frame
  min/avg/max 0.275/1.686/16.421 us; cpu 3.11%
  isr_pack        143.69/frame
    min/avg/max 6.896/8.130/10.825 us; cpu 1.87%
  isr_dma_submit  143.69/frame
    min/avg/max 0.600/0.915/1.046 us; cpu 0.21%
```

- Submission costs 0.915 µs/call versus packing 8.130 µs/call. Asynchronous wire-transfer duration is not measured by these CPU scopes.
- Inclusive `isr_wake` CPU share is 3.11%, leaving approximately 60.556 ms of foreground CPU per display window. Pack and submit run inside wake; their shares are not added. DMA-completion and other interrupts are not separately measured.
- Render already includes interrupt time. Both presets already fit one display window; no second ISR subtraction is applied.

## Summary ranking

1. `df_fused_scan` — 41.77% of frame, 26.031 ms/frame in the worst timeline window.
2. `df_lut_bake` — 25.12% of frame, 15.652 ms/frame in the worst timeline window.
3. `df_octave_noise` — 16.95% of frame, 10.564 ms/frame in the worst timeline window.
4. `df_chunk_cull` — 1.34% of frame, 0.834 ms/frame in the worst timeline window.

These costs are inclusive; octave noise is inside LUT baking. No matched native/WASM timing ledger is retained. Historical global-O3 measurements remain a separate source/configuration.

Fresh interleaved 40 s A/B captures passed both phase gates. Base: `153be2802a81baead446ff23261e1e23eab50ad9`; candidate: `182eb03a5c86145197a3723386ab0f0c66f12a3e`. Each phase used A B A B A B order, three captures per image. Variation is the larger within-image peak range; the candidate median must not exceed the base median plus that variation, with no added spills.

| Phase | Base median live peak ms | Candidate median live peak ms | Measured variation ms | Base/candidate spills |
|---|---:|---:|---:|---:|
| noise | 45.154 | 45.219 | 0.162 | 0/0 |
| balls | 37.049 | 36.817 | 0.103 | 0/0 |

Setup draws, excluded from live peaks and spills: noise base 50.171/50.108/50.175 ms; candidate 50.248/50.211/50.201 ms. balls base 19.037/19.026/19.051 ms; candidate 19.055/19.062/19.081 ms.

Qualification captures and setup/live separation: `C:/work/temp/df-core-device/final/comparison.json`. Phantasm size/layout passed; live-parent comparison found RAM1 code +176 B, unchanged RAM1 variables, FLASH code +320 B/data +60 B, no unrelated ITCM demotions, and no changed `Plot::rasterize` symbol groups.

README cells: peak 🟢 45.164 (2), spilled 🟢 0/3015 (0.00%).

README setup frame 1: 18.393 ms. Captured 2026-11-12 11:23 America/Los_Angeles.

## Caveats

- CYCCNT free-runs, so all cycle scopes include interrupt time.
- `filter_blend` first-caller parenting can hide its subtree when enabled; it is disabled here, so per-pixel scope overhead is avoided and no blend counts are claimed.
- Selective-O3 covers the ball-field evaluators and distorted-ring scan kernels under the shipping `-Os` configuration.
- No dwell, transition-speed or ordered-cycle compression was used. The explicit 300 s epoch exceeds this 190 s capture.
- Tracked and untracked source status was clean during build/capture; the artifacts retain the empty status record.
- Parser validation passed preset coverage/wrap, monotonic frames, complete counter windows, render-versus-wall and cycle/wall exactness. Mixed-parent/duplicate-name diagnostics prevent exclusive-cost interpretation.
- Qualification began with baseline A1 on COM3; after its USB failure the remaining captures and this full-cycle report used COM4, as the extraction spec permits.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=16`. `just profile DisplacementField` builds, flashes and captures through the shared-device lock; use the Setup command for this duration, window and epoch.
