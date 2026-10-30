# Millefiori on-device profile — Teensy 4.0, segmented mode (2026-10-29, **selective -O3**)

Point-in-time snapshot; this is the first canonical Millefiori capture.
Raw capture: `build/prof/millefiori_ship.log`, captured 2026-10-29 23:40 America/Los_Angeles on COM3 in the shared `C:/work/Holosphere` checkout's build directory.
Raw-log SHA-256: `68cd290451cf85f507b3d9ea612f07a1bf992e058528324bec94a05549180b31`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, live flywheel and DMA ISRs |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` hot paths, GCC 15.2.1; temporal column sampling, default effect parameters |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Millefiori, 288×144, single-entry playlist; base `9183a14385490a310dcdca16c059422d5d6482f1` plus the captured uncommitted promotion patch |
| Method | `HS_PROFILE` cycle scopes, 70 s capture, 32-frame counter windows; 11 complete windows, frames 1–352; no epoch reset or preset cycle |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-millefiori-20261029 HS_TEENSY_PORT=COM3 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Millefiori profile 70 32` |

The Phantasm show entry is 60 s. This diagnostic captures 70 s of the same default effect without a fixed-preset override, dwell compression or epoch stretch. The simulator was installed before this capture. The compiled effect has no twist lens; its authored chain retains ripple displacement, dodecahedral and tetrahedral kaleidoscopes, folded gnomonic projection, vortex warp, lattice source and two active modulation stages.

Profile image size: `FLASH: code:94992, data:153060, headers:8972, free:1774592` / `RAM1: variables:313056, code:21432, padding:11336, free:178464` / `RAM2: variables:520064, free:4224`.

Full Phantasm image size: `FLASH: code:1029504, data:840996, headers:8540, free:152576` / `RAM1: variables:312736, code:184296, padding:12312, free:14944` / `RAM2: variables:520064, free:4224`. Size/layout gates pass; FLASH total is 1,879,040 B, leaving 20,960 B under the 1,900,000 B budget.

Against the clean base's shipping Phantasm image, RAM1 code changes by **+480 B** (183,816 to 184,296), RAM1 variables by **0 B** (312,736), FLASH data by **−19,804 B** (860,800 to 840,996), and FLASH code by **+4,856 B** (1,024,648 to 1,029,504). These are whole-promotion deltas: Millefiori replaces KaleidoscopePentBright in the device roster, and Arduino-only failure-site strings retain file/line and custom messages while omitting condition text. They do not isolate Millefiori's own size.

Profile ELF SHA-256: `49d0a91a7a7e8858cc5e6705c9ed2bda8ac38bcdca9e470e9befaf2e809dc9f5`.
Phantasm ELF SHA-256: `04087e65843c8eb231e9fa2f61327255bd799f572271d2617301dfa4a0f2c78d`.
Provenance is in `build/prof/millefiori_ship.provenance`; ELF/map images, compiler/ABI attestation inputs, envdumps, build logs, `source_status.txt` and `source.diff` are retained under `build/prof/artifacts/millefiori_ship_9183a1438549_49d0a91a7a7e/`. The diff includes staged new code. The base commit alone is not the measured source. The baseline ELF and clean-source metadata are retained under `.git/fw-cache/3a072d8d5699e8b644ef918346215aef/`.

Exactness cross-check: frames **321–352** carry **3,596,238,467 cycles**, or **5,993,730.778 us**, against wall sum **5,993,731 us**: **0.037 ppm**. Every complete window is within 0.585 ppm, below the parser's 5 ppm threshold. Strict parser validation is **VALID**.

## Frame cadence

**Steady render, frame 1 excluded:** peak **116.378 ms** (frame 197), mean **101.717 ms**, **351/351 spilled frames (100.00%)**. These figures cover frames 2–352 in complete counter windows and supply the ranked README cells.

**Setup frame 1:** render **110.264 ms**, 1/1 spilled. Across all complete windows, including setup, mean render is **101.741 ms**, peak **116.378 ms**, and **352/352 frames spill**. The `fx_shader_draw` counter averages **99.590 ms/frame**; its worst window averages **108.618 ms/frame** (321–352).

The capture contains **366 complete raw frame rows**, frames 1–366. The final 14 rows (353–366) have no complete counter window and are excluded from ranked and counter aggregates. All raw rows give mean render **101.836 ms**, the same **116.378 ms** peak, and **366/366 spills**. A final partial `f` line is not a frame record.

One display window is 62.5 ms at 480 RPM. After setup, mean observed frame wall time is **187.161 ms**, about **5.34 fps** and three display windows per delivered frame. Rendering itself exceeds one display window on every measured frame; synchronization adds further idle time. `canvas_buffer_wait` is measured display-buffer idle, excluded from render telemetry. Live rendering uses the shipping quadrant clip, about **10,368 geometric pixels**, rather than the entire 41,472-pixel canvas.

## Phase-by-phase readout

This is a continuously animated, non-cycling effect. The windows below distinguish startup and observed animation cost; they are not discrete presets. Counter scopes are inclusive and nested rows overlap.

### Startup-containing window (frames 1–32)

```text
frame                    182.023 ms 109.214 Mcyc 100.00%
  pov_preserve_half       135.69 us   81.42 kcyc  0.07%
  fx_shader_draw          92.315 ms  55.389 Mcyc 50.72%
    fx_prepare_column      1.262 ms  757.02 kcyc  0.69% x144 8.76 us/col
  fx_prepare_frame        522.05 us  313.23 kcyc  0.29%
  fx_advance               1.060 ms  635.80 kcyc  0.58%
    palette_cycle_step    947.51 us  568.51 kcyc  0.52%
  fx_timeline_step         90.60 us   54.36 kcyc  0.05%
  canvas_clear             90.38 us   54.23 kcyc  0.05%
  canvas_buffer_wait      87.579 ms  52.547 Mcyc 48.11%
```

Wall min/avg/max = **88.473/182.023/198.845 ms**. This window includes setup frame 1 and the first 31 live frames. Its shader mean is lower than later high-cost portions of the animation; it is not a separate steady preset.

### Peak-frame animation window (frames 193–224)

```text
frame                    187.512 ms 112.507 Mcyc 100.00%
  pov_preserve_half       137.41 us   82.44 kcyc  0.07%
  fx_shader_draw         108.473 ms  65.084 Mcyc 57.85%
    fx_prepare_column      2.162 ms   1.297 Mcyc  1.15% x144 15.02 us/col
  fx_prepare_frame        660.05 us  396.03 kcyc  0.35%
  fx_advance              949.32 us  569.59 kcyc  0.51%
    palette_cycle_step    899.27 us  539.56 kcyc  0.48%
  fx_timeline_step        108.38 us   65.03 kcyc  0.06%
  canvas_clear             84.56 us   50.74 kcyc  0.05%
  canvas_buffer_wait      76.874 ms  46.125 Mcyc 41.00%
```

Wall min/avg/max = **180.316/187.512/196.181 ms**. Frame 197 supplies the capture's 116.378 ms render peak. Shader work dominates, while buffer idle still contributes about 76.9 ms per frame.

### Worst shader-mean window (frames 321–352)

```text
frame                    187.304 ms 112.382 Mcyc 100.00%
  pov_preserve_half       139.85 us   83.91 kcyc  0.07%
  fx_shader_draw         108.618 ms  65.171 Mcyc 57.99%
    fx_prepare_column      2.380 ms   1.428 Mcyc  1.27% x144 16.53 us/col
  fx_prepare_frame        559.39 us  335.64 kcyc  0.30%
  fx_advance              970.76 us  582.46 kcyc  0.52%
    palette_cycle_step    920.94 us  552.56 kcyc  0.49%
  fx_timeline_step        112.54 us   67.52 kcyc  0.06%
  canvas_clear             85.15 us   51.09 kcyc  0.05%
  canvas_buffer_wait      76.601 ms  45.960 Mcyc 40.90%
```

Wall min/avg/max = **181.446/187.304/193.535 ms**. Shader mean reaches 108.618 ms, slightly above the peak-frame window's 108.473 ms. Its peak individual render is 115.067 ms, showing why a worst-window mean and worst-frame peak are separate figures.

### Per-pixel figures

The direct shader assigns finished colors and does not call `filter_blend`; measured blended-pixel counts and cycles per blend are unavailable. For the worst shader-mean window, **65.171 Mcycles/frame / 10,368 geometric quadrant pixels = 6,285.75 cycles/pixel**. This includes column preparation and interrupt time; it is not an isolated per-pixel benchmark. `fx_prepare_column` runs **144 times/frame**, averaging **16.527 us/column** in that window. Temporal sampling prepares one column-time context; there are no per-pixel profiling scopes or enabled probe counts.

## Column-ISR / DMA marshaling cost

The following readout accompanies frames 321–352. Its separately timestamped **ISR window is 5,997,017 us**, including time between ISR snapshots; shares do not use the frame wall sum. Rates below divide calls by that ISR interval; calls/frame divide by 32 effect frames.

```text
isr_wake       18433.2/s  3454.50/f  0.515/1.968/20.871 us  cpu 3.63%
  isr_pack      2304.0/s   431.78/f  8.638/9.423/13.411 us  cpu 2.17%
  isr_dma_submit 2304.0/s   431.78/f  0.616/0.931/2.073 us   cpu 0.21%
```

Values after the rates are per-call min/mean/max. Pack totals **130,204 us**, versus **12,895 us** of CPU submit work in this ISR interval, about a 10.1:1 ratio. Submit measures CPU enqueue cost; asynchronous SPI wire-transfer duration is not measured by these counters.

`isr_wake` is inclusive of pack and submit. The measured flywheel ISR share is **3.63%**; the nested shares must not be added. At that share, approximately **60.231 ms** of a 62.5 ms display interval remains for foreground CPU work. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupts: the observed peak needs about **1.862×** elapsed-time speedup to reach 62.5 ms, while the worst shader-mean window's render mean needs **1.771×**. Subtracting ISR time again from render would double-count its effect.

## Summary ranking

The following hot spots use frames 321–352; percentages are of that root counter window.

1. `fx_shader_draw` — **108.618 ms/frame, 57.99%**; the direct composed pipeline dominates render work.
2. `fx_prepare_column` — **2.380 ms/frame, 1.27%**, nested within shader work; 144 temporal column preparations per frame.
3. `fx_advance` — **0.971 ms/frame, 0.52%**; includes **0.921 ms** of palette-cycle work.
4. `fx_prepare_frame` — **0.559 ms/frame, 0.30%**.

`canvas_buffer_wait` consumes 76.601 ms/frame (40.90% of root) as display idle; it is not a render hot spot. No matched native/WASM timing or global-O3 capture is available for this effect.

README cells: peak 🔴 116.378, spilled 🔴 351/351 (100.00%).

README setup frame 1: 110.264 ms. Ranked peak/spill figures exclude only setup and trailing rows without complete counter windows.

## Caveats

- Cycle scopes include ISR time and instrumentation overhead. No per-pixel scope was added.
- There is no `filter_blend` counter in this direct shader. The usual first-parent/subtree visibility artifact cannot supply missing blend counts here.
- This is the shipping selective-O3 configuration, not a global-O3 image. The path crosses `ComposedEffect::shade` and selective `Scan::Shader` hot loops; the authored placement group is `INLINE_ONLY`.
- No fixed-preset override, dwell compression or epoch stretch is used. Seventy seconds samples the continuous animation, not all possible future states.
- Captured code was uncommitted over base `9183a14385490a310dcdca16c059422d5d6482f1`; retained source status and binary diff identify the actual promotion state, including twist removal and the device roster replacement. An untracked screenshot in source status does not affect firmware.
- Setup is retained in pass aggregates and excluded from README steady figures; the entire startup-containing window remains in counter summaries. Raw trailing rows are explicitly reported separately.

Parser validation (`python tools/parse_profile.py ... validate --scope fx_shader_draw`):

```text
effect=Millefiori  windows=11
  [PASS] captured >=3 windows (11)
  [PASS] every window names one effect and resolution (Millefiori 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 11 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 11 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.0 ppm, frames 321-352)

VALID: Millefiori
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Millefiori`, `HS_PROFILE_WINDOW=32`; `just profile Millefiori` routes build, flash and capture through the shared-board lock. The exact measured 70 s command and source tree are listed above.
