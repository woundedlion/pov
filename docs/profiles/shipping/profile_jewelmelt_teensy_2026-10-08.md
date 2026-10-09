# JewelMelt on-device profile — Teensy 4.0, segmented mode (2026-10-08, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/jewelmelt_softness05_ship.log`, captured 2026-10-08 22:56 PDT on COM3.
First shipping profile for this effect. The retained capture and ELF provenance are under C:/work/temp/jewelmelt-profile-20261008.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base, newlib-nano, shipping `HS_O3` shader scan, composed shading, curl displacement and color kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | JewelMelt 288×144, softness 0.5, single-entry playlist, clean source at `22755b43d` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 70 s capture; no cadence or preset overrides |
| Reproduce | `bash tools/profile_one.sh JewelMelt profile 70 32` with `HS_PROFILE_TREE` set to the intended checkout |

Image size (`profile` env, this effect only): `FLASH: code:70108, data:155156, headers:8204` / `RAM1: variables:315040, code:19608, padding:13160, free:176480` / `RAM2: variables:520064, free:4224`.

Full Phantasm admission build: RAM1 code 177,736 B (+112 B), RAM1 variables 314,784 B (+0 B), FLASH data 847,812 B (+408 B), and FLASH code 648,404 B (+6,736 B), relative to the roster before admission. RAM1 remains 511,392 B including bank padding, leaving 12,896 B for stack; RAM2 leaves 4,224 B. All region budgets and layout invariants pass.

Exactness cross-check: frames 865–896, root counter 1,198,625,263 cycles ÷ 600 equals 1997708.772 us against 1,997,713 us measured wall sum: **2.1 ppm**. `tools/parse_profile.py ... validate` reports VALID; every window has complete frame telemetry and no epoch reset.

## Frame cadence

**Display-frame aggregate** (raw per-frame rows 2–1088): peak render is **46.39 ms** at frame 866, with **0/1087 frames (0.00%) spilled**. Setup frame 1 is excluded from both figures.

The driver renders setup frame 1 before publishing the effect to the display ISR. It sets the horizontal clip to [0,288), then installs the half-width clip and preservation hooks for subsequent draws (`hardware/pov_segmented.h`). That 77.60 ms draw initializes the full-width 288×72 segment band under the construction budget; it is not a missed 62.5 ms display deadline. `Effect::clear_buffers()` clears both buffers during construction, while `Canvas` renders and queues one buffer. The first window has 31 `pov_preserve_half` calls for 32 draws, confirming that the setup draw precedes the preservation hook.

The 62.5 ms display window gives a steady 16 fps cadence and 16.11 ms peak render margin. Live frames render one 144×72 quadrant (10,368 pixels). `canvas_buffer_wait` is the round-up idle until the next display flip. Across the raw capture, `fx_shader_draw` averages 37.62 ms/f; its worst window is 39.76 ms/f (frames 865–896).

## Phase-by-phase readout

Phase schedule: startup followed by one steady animated preset.

### Startup (window frames 1–32)

```
frame                     62.31 ms  37.39 Mcyc 100.0%
  pov_preserve_half       133.2 us   79.9 kcyc   0.2% x1.0 137.5 us/call
  fx_shader_draw          37.75 ms  22.65 Mcyc  60.6% x1.0 37748.4 us/call
  fx_prepare_frame         3.90 ms   2.34 Mcyc   6.3% x1.0 3899.8 us/call
  fx_advance               2.07 ms   1.24 Mcyc   3.3% x1.0 2070.4 us/call
  fx_timeline_step         61.9 us   37.1 kcyc   0.1% x1.0 61.9 us/call
  canvas_clear             89.7 us   53.8 kcyc   0.1% x1.0 89.7 us/call
  canvas_buffer_wait      18.28 ms  10.97 Mcyc  29.3% x1.0 18277.4 us/call
```

Wall min/avg/max = 44.56/62.31/77.60 ms. Setup frame 1 rendered 77.60 ms and remains in the raw capture and this counter tree, but is excluded from display peak and spill statistics. The window mixes that setup draw with steady rendering, so its average does not represent initialization latency.

### Steady peak (window frames 865–896)

```
frame                     62.43 ms  37.46 Mcyc 100.0%
  pov_preserve_half       136.3 us   81.8 kcyc   0.2% x1.0 136.3 us/call
  fx_shader_draw          39.76 ms  23.86 Mcyc  63.7% x1.0 39761.5 us/call
  fx_prepare_frame         3.82 ms   2.29 Mcyc   6.1% x1.0 3818.4 us/call
  fx_advance               2.06 ms   1.24 Mcyc   3.3% x1.0 2062.4 us/call
  fx_timeline_step         72.7 us   43.6 kcyc   0.1% x1.0 72.7 us/call
  canvas_clear             85.3 us   51.2 kcyc   0.1% x1.0 85.3 us/call
  canvas_buffer_wait      16.46 ms   9.88 Mcyc  26.4% x1.0 16458.7 us/call
```

Wall min/avg/max = 61.36/62.43/63.38 ms. Every frame in this window renders within 62.5 ms. Shader evaluation dominates; buffer wait fills the remaining display interval. Times and cycles are per-frame window averages; `xN` is calls per frame, with one shader draw call costing 39761.5 us.

### Per-pixel figures

This capture has no `filter_blend` counter, so blended pixels and cycles per blend are unavailable. The shader draw averages 2301.0 cycles per quadrant pixel in the steady peak window; this divides scan cost by 10,368 pixels, not by measured blend calls.

## Column-ISR / DMA marshaling cost

Steady peak window; columns are rate, min/avg/max per call, and CPU share. Shares use the separately logged ISR window interval, including report output.

```
isr_wake          1152.0/frame  0.56/1.68/19.60 us  3.08%
  isr_pack         144.0/frame  6.24/6.86/9.72 us  1.58%
  isr_dma_submit   144.0/frame  0.65/0.94/7.66 us  0.21%
```

- DMA submission averages 0.94 us/call; packing averages 6.86 us/call.
- Wire transfer runs asynchronously at the segmented driver's 24 MHz SPI clock; this capture measures submission CPU work, not transfer completion latency.
- Inclusive `isr_wake` is 3.08% CPU, leaving about 60.57 ms of foreground CPU per 62.5 ms display window. Pack and submit are nested within wake and are not added to it. Steady rendering needs no speedup to hold 16 fps. DMA-completion and other interrupts are not measured here.

## Summary ranking

1. `fx_shader_draw` — 63.7% of the steady peak window, 39.76 ms/f.
2. `fx_prepare_frame` — 6.1% of the steady peak window, 3.82 ms/f.
3. `fx_advance` — 3.3% of the steady peak window, 2.06 ms/f.
4. `pov_preserve_half` — 0.2% of the steady peak window, 0.14 ms/f.

README cells: peak 🟢 46.39, spilled 🟢 0/1087 (0.00%).

No same-workload WASM/native timing measurement accompanies this capture.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- Setup frame 1 renders the full-width segment band before display publication; display peak and spill figures use frames 2–1088, matching the archive's setup exclusion.
- `filter_blend` can appear beneath its first caller and disappear with an inactive parent; this capture emits no such counter.
- Scopes are per-frame, not per-pixel. Per-pixel scopes would add overhead and are not enabled.
- Selective-O3 is active in shader scan, composed shading, curl displacement and generated color kernels; this is the shipping `-Os` configuration.
- No dwell compression, ordered-cycle flag or epoch stretch was used. The capture stayed within the default 120 s epoch.
- The captured source tree was clean at `22755b43d`; playlist admission follows this measurement. Compiler, flags and ELF hashes are retained with the raw capture.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=JewelMelt` and `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh JewelMelt profile 70 32` builds, flashes and captures under the shared device lock.
