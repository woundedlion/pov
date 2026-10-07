# GnomonicStars on-device profile — Teensy 4.0, segmented mode (2026-10-07, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/gnomonicstars_ship.log`, captured 2026-10-07 13:42 on COM3.
Replaces `profile_gnomonicstars_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel and DMA ISRs live, COM3 |
| Image | `profile`: shipping `-Os` base and selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GnomonicStars 288×144; default Warp Speed 0.035; deterministic effect seed |
| Method | 70 s, 32-frame windows; repeat baseline in the A–B–A clock experiment |
| Source | `17390a7a1eef`: `e7f19b0d1` plus one profiling scope; original float clock |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-finding4-profile-20261007-baseline HS_TEENSY_PORT=COM3 bash tools/profile_one.sh GnomonicStars profile 70 32` |

Image size: `FLASH: code:54980, data:153136, headers:8968   free for files:1814532` / `RAM1: variables:315104, code:25352, padding:7416   free for local variables:176416` / `RAM2: variables:520064  free for malloc/new:4224`.

Exactness cross-check: frames 193–224, root cycles / 600 MHz versus wall sum differ by **1.89 ppm**; all 34 windows validate. No epoch reset or missing frame rows.

This is the unchanged-clock baseline. The [bounded-clock comparison](profile_gnomonicstars_clock_experiment_teensy_2026-10-07.md) measures an experimental implementation; it is not adopted by this report.

## Frame cadence

`gn_draw_stars` averages 9.991 ms/frame; its worst window is 15.524 ms/frame at 193–224. Peak frame render is **21.844 ms** at frame 210; **0/1087** post-setup frames spill. Average total render is 10.241608 ms. Setup frame 1 rendered 14.030 ms and is excluded from peak, spill and total-render averages. Scope averages exclude the first 32-frame window because its individual setup scope cannot be subtracted.

The display window is 62.5 ms: every captured frame holds 16 fps. The effect renders one quadrant, about 10,368 pixels. `canvas_buffer_wait` is idle until the next display flip; render is frame cost minus that wait.

## Phase-by-phase readout

The 600-star field alternates quiet footprints and larger on-screen bursts; star count stays fixed. The largest burst occurs at 193–224, with another later burst at 769–896.

### Burst (frames 193–224)

```
frame                        62.56 ms 37.54 Mcyc  100%
  pov_preserve_half          143.9 us  86.3 kcyc    0%
  gn_draw_stars              15.52 ms  9.31 Mcyc   25%
    gn_star_scan             14.82 ms  8.89 Mcyc   24% x600 25us/c
      filter_blend           739.5 us 443.7 kcyc    1% x7261 0us/c
  gn_timeline_step            39.0 us  23.4 kcyc    0%
    animation_mobius_step     10.0 us   6.0 kcyc    0%
  canvas_clear                84.4 us  50.7 kcyc    0%
  canvas_buffer_wait         46.77 ms 28.06 Mcyc   75%
```

Wall min/avg/max = 51.512/62.560/73.350 ms. Star scanning dominates render. The animation scope includes phase updates and all eight trig evaluations; its percentage is of the whole frame here.

### Quiet (frames 641–672)

```
frame                        62.54 ms 37.53 Mcyc  100%
  pov_preserve_half          145.9 us  87.6 kcyc    0%
  gn_draw_stars               9.09 ms  5.46 Mcyc   15%
    gn_star_scan              8.38 ms  5.03 Mcyc   13% x600 14us/c
      filter_blend           272.8 us 163.7 kcyc    0% x2494 0us/c
  gn_timeline_step            40.2 us  24.1 kcyc    0%
    animation_mobius_step     14.7 us   8.8 kcyc    0%
  canvas_clear                84.4 us  50.7 kcyc    0%
  canvas_buffer_wait         53.18 ms 31.91 Mcyc   85%
```

Wall min/avg/max = 59.836/62.541/65.209 ms. Smaller star footprints reduce scanning and blending; the display wait absorbs the recovered time.

### Per-pixel figures

Burst blends: 7261.0/frame, 61.1 cycles/blend. The scan consumes 1224.4 cycles per blended pixel, including non-blending scan work.

## Column-ISR / DMA marshaling cost

Burst window, per-call min/avg/max:

```
isr_wake        1154/frame  0.7/1.7/11.7 us  cpu 3.13%
isr_pack         144/frame  6.3/6.8/9.3 us  cpu 1.55%
isr_dma_submit   144/frame  0.8/0.9/1.0 us  cpu 0.21%
```

`isr_wake` is inclusive of packing and submission; those nested shares are not added. Packing is the larger CPU cost, while submission starts an asynchronous transfer. Wire-transfer duration is not measured here. At the logged wake share, nominal non-wake CPU capacity is about 60.5 ms per display window; the observed peak leaves about 40.7 ms of wall margin. DMA-completion and other interrupts are not separately measured.

## Summary ranking

1. `gn_draw_stars`: 15.524 ms/frame in the burst; the dominant render cost.
2. `pov_preserve_half`: 143.9 µs/frame.
3. `canvas_clear`: 84.4 µs/frame.
4. `gn_timeline_step`: 39.0 µs/frame; complete-capture clean-window average 43.708 µs.

No new WASM/native timing comparison was made. Native animation verification of the experimental implementation is documented in the comparison report.

README cells: peak 🟢 21.84, spilled 🟢 0/1087 (0.00%).

## Caveats

- Every cycle scope includes ISR interruptions. ISR shares use the dedicated snapshot interval.
- `filter_blend` parenting reflects its first caller; its call count is blended pixels.
- Per-star profiling has overhead common to all three experiment runs; there is no added per-pixel clock probe.
- Selective-O3 regions in the scan/filter path remain active; this is not global-O3 timing.
- Capture ran from a clean experimental branch containing one added `animation_mobius_step` scope, not a changed clock. No dwell compression or parameter override was used.
- Repeat flashing encountered a transient shared Teensy Loader busy error; the supported wrapper retry succeeded, and only the validated fresh capture is used.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=GnomonicStars`, `HS_PROFILE_WINDOW=32`. `just profile GnomonicStars` uses the same harness and the selected checkout's instrumentation.
