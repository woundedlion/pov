# LiquidWindows on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Point-in-time capture of the 60-second Phantasm slot. Raw capture `C:/work/temp/holosphere-liquid-recovery-1120/itcm.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM4 (serial 6362390), 600 MHz |
| Image | Shipping selective `-O3` regions on the `-Os` base, newlib-nano |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, synchronized, flywheel and DMA interrupts enabled |
| Effect | LiquidWindows at source `e767225183c865b8788b4b14add08c9f1d516494` |
| Method | 60 s slot, 66 s acquisition, 32 frames per counter dump |
| Reproduce | `HS_TEENSY_PORT=COM4 bash tools/profile_one.sh LiquidWindows profile 66 32` from the captured revision |

Captured 2026-11-20 17:20, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`. Raw SHA-256 `16f224b6186a017549d55a79b3c1842a90ad25dc56cf26f63d807ae8de0fd321`.

Profile image size `FLASH: code:98992, data:153068, headers:9060 free for files:1770496 / RAM1: variables:314144, code:27464, padding:5304 free for local variables:177376 / RAM2: variables:520064 free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1111928, data:724920, headers:8400 free for files:186368 / RAM1: variables:313760, code:167608, padding:29000 free for local variables:13920 / RAM2: variables:520064 free for malloc/new:4224`.

Profile ELF SHA-256 `ba1dccc7d52899888fa0d878adcf0137555f33a8b0b75790c2b5afd013584f62`. Shipping Phantasm ELF SHA-256 `ce071761435652a86282cfa5c43b31bd1840780159dc7f4e88c87df9801a16df`.

Exactness cross-check, frames 225–256. The root counter measured 1187074.552 µs, compared with 1,187,072 µs of draw time (2.15 ppm).

## Frame cadence

**Peak confirmed live render 42.228 ms** at frame 257 (LiquidWindows). **0/947 confirmed live frames spilled** over 62.5 ms. Setup frame 1 took 64.841 ms. Live render averaged 34.241 ms.

Among complete live readouts, the largest window mean was 37.160 ms/frame at frames 321–352. The dominant scope averaged 31.115 ms/frame across those readouts.

The shipping ELF assigns 480 revolutions (60 s). The single-effect harness records elapsed microseconds from its initial ISR snapshot at construction. Consecutive ISR intervals include foreground waits and USB output. Draw durations bound frame completions within the readout that crosses 60 s. Frames 2–948 are confirmed live. The raw file contains 1,032 frames, including 72 confirmed after the boundary.

Boundary resolution leaves 12 frame(s), 949–960, with uncertain membership. Their largest render was 36.408 ms and 0 exceeded the budget. The live denominator lies between 947 and 959. Ranked cells use confirmed frames. Peak and spill count hold under either cutoff.

The renderer covers one quadrant (10,368 pixels). Half a revolution is 62.5 ms, giving 16 display opportunities per second. Render telemetry subtracts `canvas_buffer_wait` from draw time. The driver also waits between draws.

README cells: peak 🟢 42.228, spilled 🟢 0/947 (0.00%).

README setup frame 1: 64.841 ms.

## Phase-by-phase readout

### Entry and setup (frames 1–32)

```
frame                    35.165 ms 21.10 Mcyc 100.0%
  pov_preserve_half      138.38 us  83.0 kcyc   0.4% x0.97 143 us/call
  fx_shader_draw         32.299 ms 19.38 Mcyc  91.9% x1 32299 us/call
    fx_prepare_column    903.98 us 542.4 kcyc   2.6% x144 6.28 us/call
  fx_palette_base_bake   812.17 us 487.3 kcyc   2.3% x1 812 us/call
    fx_palette_rgb_bake  694.24 us 416.5 kcyc   2.0% x1 694 us/call
  fx_prepare_frame       705.04 us 423.0 kcyc   2.0% x1 705 us/call
  fx_advance             928.03 us 556.8 kcyc   2.6% x1 928 us/call
    palette_cycle_step   899.87 us 539.9 kcyc   2.6% x1 900 us/call
  fx_timeline_step         9.75 us   5.8 kcyc   0.0% x1 9.75 us/call
  canvas_clear            99.94 us  60.0 kcyc   0.3% x1 99.9 us/call
  canvas_buffer_wait       5.23 us   3.1 kcyc   0.0% x1 5.23 us/call
```

Wall min/mean/max was 32.396/35.164/64.850 ms. Render averaged 35.159 ms per frame, including 32.299 ms in `fx_shader_draw`.

This regime leaves 27.341 ms of mean render headroom per display window.

### Window containing the live peak (frames 257–288)

```
frame                    36.523 ms 21.91 Mcyc 100.0%
  pov_preserve_half      141.67 us  85.0 kcyc   0.4% x1 142 us/call
  fx_shader_draw         32.639 ms 19.58 Mcyc  89.4% x1 32639 us/call
    fx_prepare_column    917.31 us 550.4 kcyc   2.5% x144 6.37 us/call
  fx_palette_base_bake    1.807 ms  1.08 Mcyc   4.9% x1 1807 us/call
    fx_palette_rgb_bake   1.690 ms  1.01 Mcyc   4.6% x1 1690 us/call
  fx_prepare_frame       806.98 us 484.2 kcyc   2.2% x1 807 us/call
  fx_advance             857.52 us 514.5 kcyc   2.3% x1 858 us/call
    palette_cycle_step   833.28 us 500.0 kcyc   2.3% x1 833 us/call
  fx_timeline_step         8.26 us   5.0 kcyc   0.0% x1 8.26 us/call
  canvas_clear            95.48 us  57.3 kcyc   0.3% x1 95.5 us/call
  canvas_buffer_wait       5.00 us   3.0 kcyc   0.0% x1 5.00 us/call
```

Wall min/mean/max was 34.301/36.523/42.232 ms. Render averaged 36.518 ms per frame, including 32.639 ms in `fx_shader_draw`.

This regime leaves 25.982 ms of mean render headroom per display window.

### Final complete in-show window (frames 897–928)

```
frame                    34.601 ms 20.76 Mcyc 100.0%
  pov_preserve_half      143.39 us  86.0 kcyc   0.4% x1 143 us/call
  fx_shader_draw         31.351 ms 18.81 Mcyc  90.6% x1 31351 us/call
    fx_prepare_column    917.58 us 550.5 kcyc   2.7% x144 6.37 us/call
  fx_palette_base_bake    1.241 ms 744.4 kcyc   3.6% x1 1241 us/call
    fx_palette_rgb_bake   1.115 ms 669.1 kcyc   3.2% x1 1115 us/call
  fx_prepare_frame       655.98 us 393.6 kcyc   1.9% x1 656 us/call
  fx_advance             936.70 us 562.0 kcyc   2.7% x1 937 us/call
    palette_cycle_step   913.45 us 548.1 kcyc   2.6% x1 913 us/call
  fx_timeline_step         8.77 us   5.3 kcyc   0.0% x1 8.77 us/call
  canvas_clear            95.52 us  57.3 kcyc   0.3% x1 95.5 us/call
  canvas_buffer_wait       5.69 us   3.4 kcyc   0.0% x1 5.69 us/call
```

Wall min/mean/max was 32.737/34.601/36.602 ms. Render averaged 34.596 ms per frame, including 31.351 ms in `fx_shader_draw`.

This regime leaves 27.904 ms of mean render headroom per display window.

### Per-pixel figures

The representative window has no complete blend counter for its 10,368 output pixels.

## Column-ISR / DMA marshaling cost

Frames 257–288, measured over the separate 1,999,618 µs interrupt interval. Each row gives calls per frame, min/mean/max cost per call, and CPU share.

```
isr_wake            1151.8/frame  0.566/1.745/19.680 us  3.21%
  isr_pack           144.0/frame  6.761/7.481/9.943 us  1.72%
  isr_dma_submit     144.0/frame  0.625/0.918/5.428 us  0.21%
```

Inclusive wake time leaves 60.494 ms per display window outside the flywheel ISR. Pack and submit run inside wake. Submission starts asynchronous SPI transfer. The 600-byte strobed column has a 230 µs transfer budget at 24 MHz. DMA completion and other interrupts have no separate timing counter.

## Summary ranking

1. `fx_shader_draw` took 32.639 ms/frame (89.4% of draw time).
2. `fx_palette_base_bake` took 1.807 ms/frame (4.9% of draw time).
3. `fx_palette_rgb_bake` took 1.690 ms/frame (4.6% of draw time).

## Regression comparison

On the Teensy, the LiquidWindows pixel shader and the veneer that calls the invariant sampler from the velocity calculation contended for the instruction cache. The recovery places the LiquidWindows pixel entry in instruction tightly coupled memory (ITCM). The timed profile image measures the full show interval. Inspection of the separately linked shipping Phantasm ELF confirms that the complete pixel entry also resides in ITCM there.

| Capture | Analyzed time | Peak live render | Spilled / confirmed live | Setup frame 1 |
|---|---:|---:|---:|---:|
| Nov 17 shipping report at `6ac5a4f04` | 30 s | 42.732 ms | 0/479 | 63.633 ms |
| Historical-source replay at `6ac5a4f04`, configured slot | 30 s | 42.723 ms | 0/472 | 63.597 ms |
| Same historical-source replay, matched duration | 60 s | 42.723 ms | 0/952 | 63.597 ms |
| Nov 20 shipping sweep at `6cca5fb7f` | 60 s | 125.926 ms | 533/533 | 213.213 ms |
| Matched pre-fix at `ce0267264` | 60 s | 125.437 ms | 533/533 | 213.132 ms |
| Recovery at `e76722518` | 60 s | 42.228 ms | 0/947 | 64.841 ms |

The pre-fix capture had 1 additional frame at the slot boundary, measuring 99.433 ms.

## Caveats

Foreground scopes include interrupt time and instrumentation overhead. The first caller determines where `filter_blend` appears. Mixed-parent and duplicate-name scope totals can overlap. The image uses shipping selective-O3 and authored preset timing. The single-entry harness measures duration from effect construction. Inter-effect handoff timing requires a playlist capture.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=LiquidWindows` and `HS_PROFILE_WINDOW=32`. Both the measured checkout and recorded source diff were clean.
