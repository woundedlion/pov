# JewelMelt on-device profile — Teensy 4.0, segmented mode (2026-11-03, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/jewelmelt_ship.log`, captured 2026-11-03 22:26 PST on COM3.
Replaces `profile_jewelmelt_selective_o3_teensy_2026-10-14.md`, which measured the preset before its retune.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os`, newlib-nano; shipping `HS_O3` shader scan, composed shading, surface displacement and color kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | JewelMelt 288×144, retuned preset (cell scale 3.0417, softness 0.2997, projection wander 0.1, central meridian 1.2566, modulation amount 0.552, pattern scale 0.6304, pattern speed 5e-5), single-entry playlist, source `2ae68ad0cdb3ee5077f3b17267786969fc6ba17a` |
| Method | `HS_PROFILE`, 32-frame windows, 70 s capture; no cadence overrides |
| Reproduce | `bash tools/profile_one.sh JewelMelt profile 70 32` with `HS_PROFILE_TREE` set explicitly |

Image size: `FLASH: code:101464, data:151632, headers:9048` / `RAM1: variables:314080, code:30072, padding:2696, free:177440` / `RAM2: variables:520064, free:4224`.

Full Phantasm passed memory/layout gates (`[teensy-gate] phantasm: PASS`): FLASH code 1,079,920 B, FLASH data 788,124 B, RAM1 code 194,472 B, RAM1 variables 313,760 B, padding 2,136 B and locals free 13,920 B.

Exactness cross-check: frames 321–352, root 3,597,700,865 cycles ÷ 600 = 5,996,168.1 us versus wall sum 5,996,167 us: **0.2 ppm**. Parser validation passed all 11 complete windows, with no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows`): `fx_shader_draw` averages 65.89 ms/f, worst window 69.14 ms/f (321–352); peak render 86.88 ms, spilled 352/352 frames, including initialization.

**Display-frame aggregate** (parser rows 2–352): peak **77.642 ms** at frame 319; **351/351 frames (100.00%) spilled**. The fastest display render is 65.022 ms, so no frame holds the 62.5 ms window. Wall time alternates between about 182 and 193 ms, so the effect displays at 5.3 fps (three display windows per frame). Fourteen trailing rows outside complete counter windows are excluded.

Setup frame 1 renders the full-width 288×72 segment band before display publication, taking 86.877 ms. Live frames render one 144×72 quadrant (10,368 pixels). The first window has 31 preservation calls for 32 draws. Setup is retained in the pass aggregate and excluded from ranked display figures. `canvas_buffer_wait` is round-up idle until the next display flip.

## Phase-by-phase readout

Phase schedule: initialization followed by one continuously animated preset.

### Startup (window frames 1–32)

```
frame                    180.45 ms 108.27 Mcyc 100.0%
  pov_preserve_half       134.4 us   80.7 kcyc   0.1% x1.0 138.8 us/call
  fx_shader_draw          63.78 ms  38.27 Mcyc  35.3% x1.0 63778 us/call
    fx_prepare_column     910.8 us  546.5 kcyc   0.5% x144 6.33 us/call
  fx_prepare_frame         3.46 ms   2.08 Mcyc   1.9% x1.0 3459.6 us/call
  fx_advance               1.14 ms  681.1 kcyc   0.6% x1.0 1135.1 us/call
    palette_cycle_step     1.04 ms  624.9 kcyc   0.6% x1.0 1041.6 us/call
  fx_timeline_step         76.6 us   46.0 kcyc   0.0% x1.0 76.6 us/call
  canvas_clear             91.4 us   54.9 kcyc   0.1% x1.0 91.4 us/call
  canvas_buffer_wait     111.54 ms  66.93 Mcyc  61.8% x1.0 111542 us/call
```

Wall min/avg/max = 68.85/180.45/189.34 ms. Frame 2 publishes without a buffer wait; every later frame waits for the third display flip. The wider setup draw sets the pass peak.

### Steady peak shader cost (window frames 321–352)

```
frame                    187.38 ms 112.43 Mcyc 100.0%
  pov_preserve_half       139.3 us   83.6 kcyc   0.1% x1.0 139.3 us/call
  fx_shader_draw          69.14 ms  41.48 Mcyc  36.9% x1.0 69140 us/call
    fx_prepare_column     947.1 us  568.3 kcyc   0.5% x144 6.58 us/call
  fx_prepare_frame         3.49 ms   2.09 Mcyc   1.9% x1.0 3488.8 us/call
  fx_advance               1.00 ms  600.8 kcyc   0.5% x1.0 1001.3 us/call
    palette_cycle_step    971.4 us  582.8 kcyc   0.5% x1.0 971.4 us/call
  fx_timeline_step         75.0 us   45.0 kcyc   0.0% x1.0 75.0 us/call
  canvas_clear             89.6 us   53.8 kcyc   0.0% x1.0 89.6 us/call
  canvas_buffer_wait     113.23 ms  67.94 Mcyc  60.4% x1.0 113229 us/call
```

Wall min/avg/max = 181.65/187.38/192.88 ms. Window shader cost ranges from 62.69 to 69.14 ms/f across the capture. Against the 2026-10-14 capture of the previous preset, `fx_shader_draw` rose from 39.61 to 69.14 ms/f while `fx_prepare_frame` and `fx_advance` fell slightly. `fx_prepare_column` is new: 144 per-column preparations per frame, about 0.95 ms/f.

### Per-pixel figures

No `filter_blend` counter is emitted, so blend counts and cycles per blend are unavailable. The steady peak shader window costs 4000.9 cycles per quadrant pixel, using 10,368 pixels rather than measured blend calls, against 2292.0 for the previous preset.

## Column-ISR / DMA marshaling cost

Steady peak shader window; shares use the separately logged 5,999,248 us ISR interval, including report output.

```
isr_wake          3455.8/frame  0.60/1.98/26.12 us  cpu 3.64%
  isr_pack         431.9/frame  8.64/9.33/14.08 us  cpu 2.14%
  isr_dma_submit   431.9/frame  0.63/0.94/3.18 us   cpu 0.21%
```

- Submission costs 0.94 us/call versus 9.33 us/call for packing. Per-frame call counts are three times the 62.5 ms case because each rendered frame spans three display windows.
- At 24 MHz SPI, the 600-byte image-plus-black composite for 72 LEDs takes 230 us including LPSPI byte framing (200 us ideal payload) asynchronously; submission timing does not measure completion latency.
- Inclusive `isr_wake` consumes 3.64% CPU, leaving about 60.23 ms foreground CPU per 62.5 ms window. Pack and submit are nested inside wake, not added to it. The 77.64 ms steady peak needs a 1.24× render speedup to hold 16 fps; the 74.15 ms steady-window mean needs 1.19×. DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 36.9% of the steady peak window, 69.14 ms/f.
2. `fx_prepare_frame` — 1.9%, 3.49 ms/f.
3. `fx_advance` — 0.5%, 1.00 ms/f, almost entirely `palette_cycle_step`.
4. `fx_prepare_column` — 0.5%, 0.95 ms/f, nested in the shader draw.
5. `pov_preserve_half` — 0.1%, 0.14 ms/f.

The previous preset peaked at 46.35 ms with zero display spills (2026-10-14, COM4, `fe0b439cb`). The retune changed seven preset values at once, and this capture does not isolate which of them raises shader cost. No same-workload WASM/native measurement accompanies this capture.

README cells: peak 🔴 77.642, spilled 🔴 351/351 (100.00%).

## Caveats

- All scopes include ISR time because CYCCNT free-runs.
- `filter_blend` normally parents under its first caller and can be hidden when that parent is inactive; this capture emits no such counter.
- Scopes are per-frame or per-column; per-pixel scope overhead is not present.
- Shipping selective-O3 covers shader scan, composed shading, surface displacement and color evaluation; this is not global-O3.
- No dwell compression, ordered cycle or epoch stretch was used; capture stayed within the default epoch.
- The main checkout carried unrelated uncommitted `docs/specs` edits; no firmware source differed from the recorded revision.
- Compiler, flags, source status and ELF/map artifacts are retained under `build/prof/artifacts/jewelmelt_ship_2ae68ad0cdb3_a3e0336bd0c5/`. Profile ELF SHA-256: `a3e0336bd0c5ef2c606ce90bfd7e1c2b91f1d322aabfa8687ade32ed186eb88d`.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=JewelMelt`, `HS_PROFILE_WINDOW=32`; `just profile JewelMelt 70` builds, flashes and captures under the shared device lock.
