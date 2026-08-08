# ChaoticStrings on-device profile — Teensy 4.0, segmented mode (2026-08-02, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ChaoticStrings`).
Raw capture: `build/prof/chaoticstrings_ship.log`. This replaces the 2026-07-25
report and is paired with
`docs/profiles/O3/profile_chaoticstrings_teensy_2026-08-02.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: GCC 15.2.1, `-Os`, newlib-nano, DMA LEDs, with selective O3 on the Plot rasterizer, geodesic strategy, and terminal blend sink |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ChaoticStrings 288×144, single-entry playlist, source snapshot `e2568f7c` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh ChaoticStrings profile 70 32` |

Image size from the exact flashed ELF: `FLASH: code:68140, data:1065232,
headers:8380` / `RAM1: variables:314496, code:49160, padding:16376,
free:144256` / `RAM2: variables:520064, free:4224`.

Compiler provenance:
`GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`.
The profile ELF SHA-256 is
`4CB41359B099E8DC07BAAF157A30BCD6C7DD95F7ABFE82894E01E4AF2F65C54A`;
its map SHA-256 is
`96D292C1C7B2EF6B7EA173B9AC1C6BF429CF63701D63586D6E8ECB69CE53C98D`.

Exactness cross-check: window frames 833–864 root counter cycles divided by
600 MHz match measured wall time within **2.1 ppm**
(`tools/parse_profile.py build/prof/chaoticstrings_ship.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `cs_build_vertices`
averages 12.628 ms/frame and whole-frame render averages 21.642 ms. The worst
window mean is 23.519 ms (frames 865–896), peak frame render is **24.854 ms**,
and **0/1088 frames (0%)** spill.

A display window is 62.5 ms; the effect renders one approximately 10,368-pixel
quadrant. The fill and steady regimes both hold 16 fps. The worst frame retains
37.646 ms of render margin. `canvas_buffer_wait` is display-flip round-up idle.

## Phase-by-phase readout

The trail grows for 115 frames, increasing vertex and raster work, then runs at
fixed history capacity. Window 97–128 covers fill completion; the peak window
is fully saturated.

### Fill completion (window frames 97–128)

```text
frame                  62.55 ms  37.53 Mcyc  100%
  cs_multiline_draw     8.84 ms   5.30 Mcyc   14%
    filter_blend       645.8 us  387.5 kcyc    1%  x6394  60.6 cyc/blend
  cs_build_vertices    12.70 ms   7.62 Mcyc   20%
  cs_noise_prepare       0.3 us     169 cyc     0%
  cs_timeline_step     118.6 us    71.2 kcyc    0%
  canvas_clear          85.3 us    51.2 kcyc    0%
  canvas_buffer_wait   40.81 ms  24.49 Mcyc   65%
```

Wall min/avg/max = 60.559/62.550/64.525 ms. Render averages 21.742 ms in
this window. Vertex construction reaches 12.70 ms as the final history frames
arrive; multiline draw follows the growing sample count at 8.84 ms.

### Saturated trail (window frames 865–896, peak of the capture)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  cs_multiline_draw    10.00 ms   6.00 Mcyc   16%
    filter_blend       760.0 us  456.0 kcyc    1%  x7594  60.0 cyc/blend
  cs_build_vertices    13.31 ms   7.99 Mcyc   21%
  cs_noise_prepare       0.2 us     145 cyc     0%
  cs_timeline_step     118.8 us    71.3 kcyc    0%
  canvas_clear          85.4 us    51.3 kcyc    0%
  canvas_buffer_wait   38.93 ms  23.36 Mcyc   62%
```

Wall min/avg/max = 60.350/62.450/64.674 ms. Render averages 23.519 ms.
Adaptive midpoint insertion makes vertex construction the dominant stage;
multiline drawing varies with the locally distorted sample count.

### Per-pixel figures

| Regime | Blends/frame | Quadrant coverage | Cycles/blend | Multiline cycles/blend |
|---|---:|---:|---:|---:|
| Fill completion | 6,394 | 0.617× | 60.6 | 829.0 |
| Saturated peak | 7,594 | 0.732× | 60.0 | 790.3 |

`filter_blend` counts antialiased pixel writes, including repeated blends to a
pixel. The multiline figure includes edge setup, adaptive geodesic sampling,
shading, and the terminal blend.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.61/1.87/11.26 us  cpu 3.44%
isr_pack         144/frame  min/avg/max 6.11/6.61/9.01 us   cpu 1.52%
isr_dma_submit   144/frame  min/avg/max 0.83/0.94/2.71 us   cpu 0.21%
```

- Packing costs about 7× submit per call; DMA wire transfer is asynchronous
  and excluded from these CPU timings.
- The three ISR groups consume 5.17% of wall time in the peak window.
- Neither phase needs a cadence speedup; the peak is 37.646 ms below one
  62.5 ms display window.

## Summary ranking

1. `cs_build_vertices` — 13.31 ms, 56.6% of saturated render: warp and adaptive
   midpoint evaluation across the full trail.
2. `cs_multiline_draw` — 10.00 ms, 42.5%: geodesic rasterization and 7,594
   antialiased blend writes per frame.
3. `cs_timeline_step` — 0.119 ms: orientation, motion, and palette drivers.

The prior 2026-07-25 image peaked at 16.14 ms before adaptive distortion
sampling. The smoother current path adds workload but remains comfortably in
the 16 fps tier.

## Caveats

- Cycle scopes include ISR time because CYCCNT free-runs.
- `filter_blend` parents under `cs_multiline_draw`; its calls count writes, not
  unique pixels.
- Selective O3 reaches Plot rasterization, geodesic interpolation, and the
  terminal pipeline blend sink.
- Per-pixel scope overhead is included in `filter_blend`.
- The capture used committed effect code at `e2568f7c`; the source tree also
  contained pre-existing documentation edits, all outside measured code paths.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ChaoticStrings`,
`HS_PROFILE_WINDOW=32`; the reproduce command builds, locks, flashes, captures,
validates, and archives the exact artifacts.

