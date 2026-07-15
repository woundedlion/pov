# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-07-15, **-O3**)

Point-in-time snapshot at the slim RingGroup tip (`730b6d2f` — inline
`stroke_alpha` + covering-ring row intervals); supersedes the fused-only
tip (`56d8c854`) and the 2026-07-14 legacy report. Raw capture:
`build/profile_capture_ringspin_slim_o3_w32.log` (32-frame windows; port
dropped at frame 768, ~1.6 epochs — ample). The global-`-O3` ceiling twin of the shipping
[selective-O3 report](../shipping/profile_ringspin_teensy_2026-07-15.md);
setup identical except the optimization level.

**Headline: `-O3` + the slim scan locks RingSpin at a hard 16 fps — every
window, including the trail-ramp window the fused-only `-O3` could not
hold.** All 24 captured windows average 60.4–62.8 ms wall; the slim
per-blend eval (502 cyc/blend, down from the fused-only tip's 665) pulls
the peak-coverage ramp frame under the 62.5 ms window.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: same Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE -D HS_PROFILE_WINDOW=32`, but base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `RingSpin<288, 144>` at master `730b6d2f` — slim `Scan::RingGroup` group scan (76 passes/frame) |
| Method | DWT `CYCCNT` scopes + `micros()` wall per `draw_frame`, dumped every 32 frames |
| Reproduce | `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 80,028 B; ITCM (RAM1 code) 64,392 B; RAM2 free
4,736 B. (At -Os: 42,228 / 28,792 / 4,736 — FLASH +89%, ITCM +124%.)

**Exactness cross-check** — window frames 97–128: root counter
1,200,400,418 cyc = 2,000,667.4 µs vs measured wall sum 2,000,671 µs
(Δ ≈ 1.8 ppm).

## Frame cadence

Same 62.5 ms display-window quantization as the `-Os` twin (quadrant ≈
10,368 px per `draw_frame`). **Every captured window locks 16 fps**: wall
avg 60.4–62.8 ms across all 24 windows, worst single frame 101 ms (the
epoch-start filling frame). Render spans 31–42 ms — always under the window,
including the trail-ramp window (frames 33–64: render 41.8 ms, wall 61.6 ms),
which the fused-only `-O3` ran at 118.8 ms.

## Regime-by-regime readout (32-frame windows, per-frame averages)

### Heavy coverage (frames 97–128)

```
frame                   62.52 ms  37.51 Mcyc  100%
  rs_draw_rings         39.42 ms  23.65 Mcyc   63%
    rs_ring_scan        39.06 ms  23.44 Mcyc   62%  x76     514 us/group
      filter_blend       8.01 ms   4.81 Mcyc   13%  x46681  103 cyc/blend
  rs_buffer_wait        22.83 ms  13.70 Mcyc   37%
  rs_timeline_step      50.6 us   30.4 kcyc     0%
```

Wall min 54.2 / avg 62.5 / max 72.2 ms. The heaviest sustained coverage
(46.7 k blends/frame) renders in 39.1 ms at **502 cyc per blended pixel**
(665 at the fused-only tip, **−24%**), leaving 23 ms of window slack.

### Light coverage (frames 193–224)

```
frame                   62.72 ms  37.63 Mcyc  100%
  rs_draw_rings         31.67 ms  19.00 Mcyc   50%
    rs_ring_scan        31.29 ms  18.78 Mcyc   50%  x76     412 us/group
      filter_blend       6.57 ms   3.94 Mcyc   10%  x37920  104 cyc/blend
  rs_buffer_wait        30.80 ms  18.48 Mcyc   49%
  rs_timeline_step      52.0 us   31.2 kcyc     0%
```

Wall min 54.1 / avg 62.7 / max 70.7 ms — 495 cyc/blend; render 31.3 ms.

### Trail ramp — now locked (frames 33–64)

```
frame                   61.56 ms  36.94 Mcyc  100%
  rs_draw_rings         41.77 ms  25.06 Mcyc   68%
    rs_ring_scan        41.77 ms  25.06 Mcyc   68%  x76     550 us/group
      filter_blend       8.96 ms   5.37 Mcyc   15%  x52199  103 cyc/blend
  rs_buffer_wait        19.15 ms  11.49 Mcyc   31%
  rs_timeline_step      52.4 us   31.4 kcyc     0%
```

Wall min 28.1 / avg 61.6 / max 97.8 ms. Peak trail-fill coverage
(52.2 k blends/frame) renders in 41.8 ms — **the fused-only `-O3` ran this
same window at 118.8 ms (68 ms render)**, so this is exactly the frame the
slim per-blend eval rescued into a lock.

### Per-pixel figures

**495–502 cyc (0.83 µs) per blended pixel** across regimes, down from the
fused-only tip's 665 (−24%); `filter_blend` ~103 cyc/blend (107 before).

## Column-ISR / DMA marshaling cost

All-locked windows: `isr_wake` 2.4–3.0%, `isr_pack` 1.0–1.5%,
`isr_dma_submit` 0.16–0.21% CPU. ~1.7 ms ISR per 62.5 ms window ⇒ ~60.8 ms
render budget; the heaviest render (39 ms) clears it with >20 ms margin.

## Summary ranking (heavy window, share of the 62.5 ms frame)

1. `rs_ring_scan` (slim `Scan::RingGroup::draw` × 76) — **62%** (39.1 ms)
2. display-window sync (`rs_buffer_wait`) — 37% (22.8 ms)
3. `filter_blend` (inside the scan) — 13% (8.0 ms)
4. everything else — <0.1%

## -O3 vs -Os (both at the slim tip `730b6d2f`)

| metric | -Os | -O3 | Δ |
|---|---|---|---|
| `rs_ring_scan`, heavy | 49.6 ms | 39.1 ms | 1.27× faster |
| scan cyc / blended px | 700 | 502 | 1.39× |
| `filter_blend` per blend | 149 cyc | 103 cyc | 1.45× |
| cadence | ~14 fps effective | **hard-locked 16 fps** | tier gained |
| image (FLASH / ITCM) | 42.2 / 28.8 KB | 80.0 / 64.4 KB | +89% / +124% |

Slim-vs-fused at `-O3`: scan 51.6 → 39.1 ms, 665 → 502 cyc/blend (−24%),
and the cadence goes from "16 fps with one 2-window ramp per epoch" to a
**hard lock**. Combined with the `-Os` result, RingSpin is a clean
selective-`-O3` candidate: the slim scan already lifts `-Os` to ~14 fps,
and `-O3` on just the fused-scan translation units finishes the lock.

## Caveats

- Not the shipping config: the full 26-effect roster overflows ITCM at
  `-O3`; single-effect image only.
- All scopes absorb live ISR time; `filter_blend` parents under
  `rs_ring_scan`.
- `-ffast-math` changes float semantics vs the `-Os`/host builds; output
  is not bit-comparable across optimization levels.
- Same slim covering-ring AA-tail note as the `-Os` twin.
- Capture ended at frame 768 (USB serial drop); the 24 dumped windows cover
  ~1.6 epochs and every regime, all locked.

## Harness

`pio run -e profile_o3 -t upload` with
`PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=RingSpin -D HS_PROFILE_WINDOW=32'`,
then `python tools/profile_capture.py --seconds 150 --out <log>`.
