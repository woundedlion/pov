# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-11-07, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `build/prof/displacementfield_ship.log`. Captured 2026-11-07 23:11 America/Los_Angeles on COM3 (serial 6362410). Replaces `profile_displacementfield_teensy_2026-11-05.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, newlib-nano, selective `HS_O3`; `DisplacementField::ball_field`, `ball_field_at`, and `Scan::DistortedRingStack` draw kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, source `c79efaa162384623b6eec0778df050c59142b736` |
| Method | One 150 s pass, 16-frame counter windows, normal preset timing and default 3,600 s profile epoch; no compression knobs |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere bash tools/profile_one.sh DisplacementField profile 150 16` |

Image size:

```text
teensy_size:   FLASH: code:116528, data:150508, headers:8420   free for files:1756160
teensy_size:    RAM1: variables:314464, code:25112, padding:7656   free for local variables:177056
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `da9f59eab04f962691de0cc8cd0e8e8a1472a2c2a6d30c40251f44d7ef1158b3`. Phantasm ELF SHA-256: `ccf2a298a4a9f78dbc32dd406f2e3e2209f02ca90c795dfd6a58a3a5803b6810`. Both retained artifact hashes were verified against capture provenance. The adjacent `.provenance`, envdump and build logs retain compiler/ABI and build flags; artifact directory: `build/prof/artifacts/displacementfield_ship_c79efaa16238_da9f59eab04f`.

Exactness cross-check: frames 2113–2128, root 598,095,351 cycles ÷ 600 MHz = 996825.585 µs versus measured wall sum 996,824 µs: **1.590 ppm** difference, below the 5 ppm gate.

Both `phantasm` and `profile` builds passed their size/layout gates, and the wrapper passed compiler/ABI attestation and fresh-image checks. A vendor `usb_serial.c` signedness warning occurred in the profile build. No engine source changed, so native behavior tests were not rerun.

## Frame cadence

**Live render, setup frame 1 excluded:** peak **45.107 ms** (frame 603), **0/2367 spilled** (0.00%). Setup frame 1: **18.485 ms**. Complete counter windows cover frames 1–2368; the incomplete trailing window is excluded. Live transitions remain included.

Pass aggregate (`parse_profile.py ... windows`): peak 45.11 ms, 0/2368 spilled including setup. Live mean render is 36.628 ms/frame. `df_timeline_step` averages 36.327 ms/frame; its worst window is 44.006 ms/frame (frames 2113–2128).

At 480 RPM a half revolution is 62.5 ms, yielding 16 fps when rendering fits one window. The effect renders one quadrant, 72×144 = 10,368 pixels. Both noise and balls, including their transitions, fit 16 fps: noise peak leaves 17.393 ms and balls leaves 24.957 ms. `canvas_buffer_wait` is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 150-frame noise fade-in, 600-frame hold, 150-frame fade-out; then 900 spawning frames of balls followed by their drain, returning to noise. Palette wipes and ambient motion continue. Preset markers confirm noise → balls → noise within this capture.

Counter trees preserve logged nesting and every duplicate counter instance. Tagged mixed-parent or duplicate-name scopes overlap; they must not be summed as exclusive costs. In particular, `ring_stack_table` can appear as a separate root because its first caller varies.

### Startup / noise entry (frames 1–16)

```text
frame                          57.750 ms 34.650 Mcyc 100.0%
  pov_preserve_half           137.255 us 82.353 kcyc   0.2%
    x0.94/frame; 146.405 us/call
  df_timeline_step             25.081 ms 15.049 Mcyc  43.4%
    df_draw_rings              24.991 ms 14.995 Mcyc  43.3%
      df_fused_scan             9.842 ms  5.905 Mcyc  17.0%
        x0.56/frame; 17497.669 us/call
      df_lut_bake               7.001 ms  4.201 Mcyc  12.1%
        df_hue_table_prep     202.512 us 121.507 kcyc   0.4%
          x7.94/frame; 25.513 us/call
        df_octave_noise         4.870 ms  2.922 Mcyc   8.4%
          x26.81/frame; 181.642 us/call
      df_chunk_cull           511.727 us 307.036 kcyc   0.9%
        x26.88/frame; 19.041 us/call
      df_lut_bake               1.138 ms 682.913 kcyc   2.0%
        df_hue_table_prep      10.375 us  6.225 kcyc   0.0%
          x5.25/frame; 1.976 us/call
      df_chunk_cull           268.750 us 161.250 kcyc   0.5%
        x18.00/frame; 14.931 us/call
      df_fused_scan             4.885 ms  2.931 Mcyc   8.5%
        ring_stack_table        1.576 ms 945.502 kcyc   2.7%
          x1.50/frame; 1050.558 us/call
        ring_stack_table       87.382 us 52.429 kcyc   0.2%
          x0.06/frame; 1398.110 us/call
    animation_noise_clock       0.165 us   99.0 cyc    0.0%
      x1.00/frame; 0.165 us/call
  canvas_clear                100.472 us 60.283 kcyc   0.2%
    x1.00/frame; 100.472 us/call
  canvas_buffer_wait           32.398 ms 19.439 Mcyc  56.1%
    x1.00/frame; 32397.604 us/call
  df_prepare_fields             1.111 us  666.5 cyc    0.0%
    x1.00/frame; 1.111 us/call
```

Wall min/avg/max = 15.292/57.750/79.861 ms. Window mean render is 25.353 ms/frame. The first draw costs 18.485 ms; later draws in this window include the noise fade-in. Wall startup transients include display synchronization and do not indicate a render spill.

### Held noise, worst timeline window (frames 2113–2128)

```text
ring_stack_table                2.092 ms  1.255 Mcyc   3.4%
  x2.00/frame; 1046.177 us/call
frame                          62.302 ms 37.381 Mcyc 100.0%
  pov_preserve_half           143.619 us 86.171 kcyc   0.2%
    x1.00/frame; 143.619 us/call
  df_timeline_step             44.006 ms 26.404 Mcyc  70.6%
    df_draw_rings              43.894 ms 26.336 Mcyc  70.5%
      df_fused_scan            26.051 ms 15.630 Mcyc  41.8%
        x1.00/frame; 26050.520 us/call
      df_lut_bake              15.325 ms  9.195 Mcyc  24.6%
        df_hue_table_prep       1.181 ms 708.306 kcyc   1.9%
          x23.56/frame; 50.101 us/call
        df_octave_noise        10.416 ms  6.250 Mcyc  16.7%
          x48.00/frame; 217.006 us/call
      df_chunk_cull           950.271 us 570.162 kcyc   1.5%
        x48.00/frame; 19.797 us/call
    animation_noise_clock       0.160 us   96.0 cyc    0.0%
      x1.00/frame; 0.160 us/call
  canvas_clear                 95.646 us 57.387 kcyc   0.2%
    x1.00/frame; 95.646 us/call
  canvas_buffer_wait           18.017 ms 10.810 Mcyc  28.9%
    x1.00/frame; 18016.633 us/call
  df_prepare_fields             2.289 us  1.373 kcyc   0.0%
    x1.00/frame; 2.289 us/call
```

Wall min/avg/max = 61.346/62.301/63.670 ms. Window mean render is 44.285 ms/frame. Noise LUT work includes octave sampling; fused scan remains the largest draw scope. This is the highest window mean, while the single-frame live maximum occurs in frames 593–608.

### Held balls (frames 1617–1632)

```text
ring_stack_table                1.747 ms  1.048 Mcyc   2.8%
  x2.00/frame; 873.453 us/call
frame                          62.345 ms 37.407 Mcyc 100.0%
  pov_preserve_half           142.192 us 85.315 kcyc   0.2%
    x1.00/frame; 142.192 us/call
  df_timeline_step             35.926 ms 21.556 Mcyc  57.6%
    df_draw_rings              35.695 ms 21.417 Mcyc  57.3%
      df_fused_scan            20.692 ms 12.415 Mcyc  33.2%
        x1.00/frame; 20691.949 us/call
      df_lut_bake              12.838 ms  7.703 Mcyc  20.6%
        df_hue_table_prep     369.772 us 221.863 kcyc   0.6%
          x36.75/frame; 10.062 us/call
      df_chunk_cull           655.668 us 393.401 kcyc   1.1%
        x44.62/frame; 14.693 us/call
    animation_noise_clock       0.160 us   96.0 cyc    0.0%
      x1.00/frame; 0.160 us/call
  canvas_clear                 94.438 us 56.663 kcyc   0.2%
    x1.00/frame; 94.438 us/call
  canvas_buffer_wait           26.105 ms 15.663 Mcyc  41.9%
    x1.00/frame; 26104.876 us/call
  df_prepare_fields            18.641 us 11.185 kcyc   0.0%
    x1.00/frame; 18.641 us/call
```

Wall min/avg/max = 60.972/62.345/63.691 ms. Window mean render is 36.240 ms/frame. The ball path avoids the octave-noise work of the noise phase. Its lower render cost leaves more time in the display-buffer wait.

### Noise-to-balls handoff (frames 897–912)

```text
frame                          62.616 ms 37.570 Mcyc 100.0%
  pov_preserve_half           147.802 us 88.681 kcyc   0.2%
    x1.00/frame; 147.802 us/call
  df_timeline_step             14.185 ms  8.511 Mcyc  22.7%
    df_draw_rings              14.049 ms  8.430 Mcyc  22.4%
      df_fused_scan             8.930 ms  5.358 Mcyc  14.3%
        x0.75/frame; 11907.241 us/call
      df_lut_bake             975.932 us 585.559 kcyc   1.6%
        df_hue_table_prep      24.380 us 14.628 kcyc   0.0%
          x3.88/frame; 6.292 us/call
      df_chunk_cull            85.086 us 51.052 kcyc   0.1%
        x3.88/frame; 21.958 us/call
      df_lut_bake             524.766 us 314.859 kcyc   0.8%
        df_hue_table_prep       5.454 us  3.272 kcyc   0.0%
          x3.00/frame; 1.818 us/call
      df_chunk_cull           121.496 us 72.898 kcyc   0.2%
        x8.75/frame; 13.885 us/call
      df_fused_scan             2.463 ms  1.478 Mcyc   3.9%
        ring_stack_table        1.458 ms 874.595 kcyc   2.3%
          x1.75/frame; 832.947 us/call
    animation_noise_clock       0.162 us   97.1 cyc    0.0%
      x1.00/frame; 0.162 us/call
  canvas_clear                 89.234 us 53.540 kcyc   0.1%
    x1.00/frame; 89.234 us/call
  canvas_buffer_wait           48.136 ms 28.882 Mcyc  76.9%
    x1.00/frame; 48136.183 us/call
  df_prepare_fields             9.581 us  5.749 kcyc   0.0%
    x1.00/frame; 9.581 us/call
```

Wall min/avg/max = 59.039/62.616/63.751 ms. Window mean render is 14.480 ms/frame. The fading noise and newly adopted balls produce a mixed window, excluded from clean-hold selection. It remains included in live peaks and spills, and has substantially less render work than either held maximum.

### Per-preset table

Preset 1 is `noise`; preset 2 is `balls`. Serial markers `Preset: 2/2` then `Preset: 1/2` confirm return to the initial noise preset (zero-based 0). The parser infers the initial noise ownership. Clean holds exclude mixed-owner windows and gate on the modal full-frame `df_timeline_step` call count (16 per window). Rows use the costliest qualifying window and rank by timeline cost.

| # / preset | Windows / clean | Blended px/f | Timeline ms/f | Render ms/f | fps |
|---|--:|--:|--:|--:|--:|
| 1 / noise | 85/84 | unavailable | 44.006 | 44.285 | 16.05 |
| 2 / balls | 63/62 | unavailable | 35.926 | 36.240 | 16.04 |

Cadence buckets include live transitions and exclude only setup frame 1; these are stricter than clean-hold selection.

| Preset | Peak live render ms | Spilled/frames | Cadence |
|---|--:|--:|---|
| Noise | 45.107 | 0/1355 | 16 fps |
| Balls | 37.543 | 0/1012 | 16 fps |

### Per-pixel figures

Per-pixel blend scopes are disabled in the shipping timing image. Blended pixels/frame, coverage, cycles/blend and scan cycles/blended pixel are unavailable; quadrant area is not substituted for actual blended work.

## Column-ISR / DMA marshaling cost

Last complete window, 16 frames; separately measured ISR interval 1,000,901 µs:

```text
isr_wake        1153.12/frame
  min/avg/max 0.268/1.660/25.081 us; cpu 3.06%
  isr_pack      144.12/frame
    min/avg/max 7.001/8.011/10.736 us; cpu 1.84%
  isr_dma_submit 144.12/frame
    min/avg/max 0.598/0.915/1.108 us; cpu 0.21%
```

- Submission costs 0.915 µs/call versus packing's 8.011 µs/call; asynchronous wire-transfer duration is not measured by these CPU scopes.
- Inclusive `isr_wake` CPU share is 3.06%, leaving approximately 60.588 ms of foreground CPU per 62.5 ms window. Pack and submit run inside wake; their shares are not added to it. DMA-completion and other interrupts are not separately measured.
- Render already includes interrupt time. Noise and balls each require 1.000× speedup to fit one display window (both already fit); no second ISR subtraction is applied to measured render.

## Summary ranking

1. `df_fused_scan` — 41.81% of frame, 26.051 ms/frame in the worst timeline window.
2. `df_lut_bake` — 24.60% of frame, 15.325 ms/frame in the worst timeline window.
3. `df_octave_noise` — 16.72% of frame, 10.416 ms/frame in the worst timeline window.
4. `df_chunk_cull` — 1.53% of frame, 0.950 ms/frame in the worst timeline window.

Costs are inclusive: octave noise is inside LUT baking. The new live peak is +0.196 ms against the previous 44.911 ms snapshot, with zero spills in both; source and board differ, so this is not a controlled A/B regression claim. Historical global-O3 peak is 57.45 ms; no matched native/WASM ledger is retained in this checkout.

README cells: peak 🟢 45.107 (2), spilled 🟢 0/2367 (0.00%).

README setup frame 1: 18.485 ms. Captured 2026-11-07 23:11 America/Los_Angeles.

## Caveats

- All cycle scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting can hide its subtree under an inactive parent when enabled. It is disabled here; per-pixel scope overhead is avoided and no blend counts are claimed.
- Shipping selective-O3 regions cover ball-field evaluators and the distorted ring-stack scan. This is the shipping `-Os` configuration, not a global-O3 comparison.
- No preset dwell, transition-speed or ordered-cycle override was used. The default one-hour profile epoch prevents reset during the 150 s capture.
- Tracked source was clean during build/capture. Two unrelated untracked specifications were present: `docs/specs/phantasm_wide_canvas.md` and `docs/specs/test_machinery_simplification.md`; neither participates in firmware compilation. Source status and tracked diff are retained beside the ELF artifacts.
- Parser validation passed every check across 148 windows, including monotonic frame numbering, all presets, wrap, complete telemetry, render-versus-wall and cycle/wall exactness. Mixed-parent/duplicate-name diagnostics are informational and prevent exclusive-cost interpretation.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=16`; default `HS_PROFILE_EPOCH_REVS=28800` at 480 RPM. `just profile DisplacementField` builds, flashes and captures through the shared-device lock. Use the Setup reproduction command for this full-cycle duration and 16-frame resolution.
