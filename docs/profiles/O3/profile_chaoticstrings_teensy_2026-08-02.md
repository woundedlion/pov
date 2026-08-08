# ChaoticStrings on-device profile — Teensy 4.0, segmented mode (2026-08-02, **-O3**)

Point-in-time global-O3 reference paired with the
[`shipping report`](../shipping/profile_chaoticstrings_teensy_2026-08-02.md).
Raw capture: `build/prof/chaoticstrings_o3.log`. This is the first current O3
twin for ChaoticStrings.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz on COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile_o3` env: GCC 15.2.1, global `-O3 -ffast-math -fno-finite-math-only`, newlib-nano, DMA LEDs |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ChaoticStrings 288×144, single-entry playlist, source snapshot `e2568f7c` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh ChaoticStrings profile_o3 70 32` |

Image size from the exact flashed ELF: `FLASH: code:96596, data:1065048,
headers:8780` / `RAM1: variables:314496, code:69848, padding:28456,
free:111488` / `RAM2: variables:520064, free:4224`.

Compiler provenance:
`GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`.
The profile ELF SHA-256 is
`0AA377553E204E6FF2106F61E87A7064FE95B7512E633969688F62BF0A5800A8`;
its map SHA-256 is
`8CF69A709C858B58D9B010B6739FC0772B995C47DA7CCA0DC635A987A37E8DD8`.

Exactness cross-check: window frames 161–192 root counter cycles divided by
600 MHz match measured wall time within **2.1 ppm**
(`tools/parse_profile.py build/prof/chaoticstrings_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `cs_build_vertices`
averages 11.516 ms/frame and whole-frame render averages 19.602 ms. The worst
window mean is 21.028 ms (frames 1025–1056), peak frame render is **22.157 ms**
(frames 865–896), and **0/1088 frames (0%)** spill.

The 62.5 ms display window is held throughout fill and steady state. The peak
retains 40.343 ms of margin while rendering one approximately 10,368-pixel
quadrant. `canvas_buffer_wait` is display-flip round-up idle.

## Phase-by-phase readout

The first 115 frames fill the trail; subsequent frames run at fixed history
capacity with distortion-dependent midpoint insertion.

### Fill completion (window frames 97–128)

```text
frame                  62.55 ms  37.53 Mcyc  100%
  cs_multiline_draw     7.92 ms   4.75 Mcyc   13%
    filter_blend       594.3 us  356.6 kcyc    1%  x6394  55.8 cyc/blend
  cs_build_vertices    11.58 ms   6.95 Mcyc   19%
  cs_noise_prepare       0.2 us     128 cyc     0%
  cs_timeline_step      91.4 us    54.8 kcyc    0%
  canvas_clear          85.0 us    51.0 kcyc    0%
  canvas_buffer_wait   42.87 ms  25.72 Mcyc   69%
```

Wall min/avg/max = 60.517/62.550/64.548 ms. Render averages 19.677 ms.
Global O3 reduces both dominant stages as the trail reaches capacity.

### Saturated trail (window frames 865–896, peak of the capture)

```text
frame                  62.45 ms  37.47 Mcyc  100%
  cs_multiline_draw     8.67 ms   5.20 Mcyc   14%
    filter_blend       665.8 us  399.5 kcyc    1%  x7174  55.7 cyc/blend
  cs_build_vertices    12.16 ms   7.30 Mcyc   19%
  cs_noise_prepare       0.2 us     125 cyc     0%
  cs_timeline_step      95.2 us    57.1 kcyc    0%
  canvas_clear          85.2 us    51.1 kcyc    0%
  canvas_buffer_wait   41.44 ms  24.86 Mcyc   66%
```

Wall min/avg/max = 60.451/62.449/64.376 ms. Render averages 21.009 ms.
The peak is still vertex-build dominated. Fast-math changes a few adaptive
sample decisions, so the blend count is not bit-identical to shipping.

### Per-pixel figures

| Regime | Blends/frame | Quadrant coverage | Cycles/blend | Multiline cycles/blend |
|---|---:|---:|---:|---:|
| Fill completion | 6,394 | 0.617× | 55.8 | 743.4 |
| Saturated peak | 7,174 | 0.692× | 55.7 | 725.0 |

Counts are antialiased writes rather than unique pixels. Multiline cycles
include adaptive sampling and geodesic raster setup as well as blending.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.53/1.78/11.02 us  cpu 3.27%
isr_pack         144/frame  min/avg/max 6.34/6.89/9.09 us   cpu 1.58%
isr_dma_submit   144/frame  min/avg/max 0.74/0.91/1.03 us   cpu 0.20%
```

- Packing dominates CPU-side DMA preparation; wire transfer is asynchronous
  and excluded from these scopes.
- The three ISR groups consume 5.05% of wall time in the peak window.
- No phase needs a cadence speedup; peak render is 40.343 ms below one display
  window.

## Summary ranking

1. `cs_build_vertices` — 12.16 ms, 57.9% of saturated render.
2. `cs_multiline_draw` — 8.67 ms, 41.3% of saturated render.
3. `cs_timeline_step` — 0.095 ms.

## Caveats

- Cycle scopes include live ISR execution.
- This single-effect global-O3 image is a compiler ceiling, not a shippable
  full-roster configuration.
- `filter_blend` parents under `cs_multiline_draw`; calls count writes rather
  than unique pixels, and per-pixel scope overhead is included.
- Fast-math perturbs adaptive subdivision thresholds and therefore exact
  raster workload.
- The capture used committed effect code at `e2568f7c`; unrelated documentation
  changes were present in the source tree but outside measured code paths.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ChaoticStrings`,
`HS_PROFILE_WINDOW=32`; use the reproduce command for locked build, flash,
capture, validation, and artifact archival.

## Global O3 vs selective O3

Global O3 peaks at 22.157 ms versus 24.854 ms for shipping, a **2.697 ms
(10.9%) improvement**. In the saturated peak window, vertex construction falls
13.31 → 12.16 ms and multiline drawing falls 10.00 → 8.67 ms. FLASH code grows
by **28,456 B** and ITCM code by **20,688 B**, with no cadence change.
