# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-07-15, **-O3**)

Point-in-time snapshot at the landed fused trail-scan tip (`56d8c854`);
replaces the 2026-07-14 legacy-path report. Raw capture:
`build/profile_capture_ringspin_fused_o3_w32.log` (32-frame windows,
~150 s, 2 epochs). This is the **-O3** twin of the shipping `-Os` report
([../Os/profile_ringspin_teensy_2026-07-15.md](../Os/profile_ringspin_teensy_2026-07-15.md));
setup is identical except the optimization level.

**Headline: `-O3` + the fused scan locks RingSpin at 16 fps.** Every
32-frame window of both epochs averages 62.1–62.9 ms wall except the one
trail-ramp window per epoch (frames 33–64, 118.8 ms) — effective
≈ **15.6 fps**, up from the baseline's borderline 16↔8 oscillation.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: same Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE -D HS_PROFILE_WINDOW=32`, but base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `RingSpin<288, 144>` at master `56d8c854` — each trail frame's ≤4 sub-rings drawn as one fused `Scan::RingGroup` pass (76 group passes/frame) |
| Method | DWT `CYCCNT` scopes + `micros()` wall per `draw_frame`, dumped every 32 frames |
| Reproduce | `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 86,236 B; ITCM (RAM1 code) 70,424 B; RAM2 free
4,736 B. (At -Os: 41,908 / 28,472 / 4,736 — FLASH +106%, ITCM +147%; the
fused group loop inlines aggressively at -O3.)

**Exactness cross-check** — window frames 33–64 (epoch 1): root counter
2,280,179,716 cyc = 3,800,299.5 µs vs measured wall sum 3,800,301 µs
(Δ ≈ 0.4 ppm).

## Frame cadence

Same 62.5 ms display-window quantization as the `-Os` twin (one quadrant ≈
10,368 px per `draw_frame`). Observed:

- **Everything after the trail ramp**: **62.5 ms (16 fps), locked** — wall
  avg 62.1–62.9 ms in every window of both epochs, worst single frame
  82.9 ms. Render spans 40–52 ms, always under the window.
- **Trail ramp + early hot orientations** (frames 33–64 of each epoch,
  bit-identical across epochs): render 69.0 ms → 2-window frames,
  window avg 118.8 ms.

`rs_buffer_wait` is the round-up idle; in locked windows it is the 10–21 ms
residual to the window boundary.

## Regime-by-regime readout (32-frame windows, per-frame averages)

### Locked 16 fps, heavy coverage (frames 97–128)

```
frame                   62.38 ms  37.43 Mcyc  100%
  rs_draw_rings         52.12 ms  31.27 Mcyc   84%
    rs_ring_scan        51.55 ms  30.93 Mcyc   83%  x76     678 us/group
      filter_blend       8.31 ms   4.99 Mcyc   13%  x46520  107 cyc/blend
  rs_buffer_wait        10.21 ms   6.13 Mcyc   16%
  rs_timeline_step      50.8 us   30.5 kcyc     0%
```

Wall min 52.9 / avg 62.4 / max 74.0 ms. The heaviest sustained coverage
(46.5 k blends/frame — the baseline's "in-band" regime) renders in 52 ms,
leaving 10 ms of window slack; scan cost is **665 cyc per blended pixel**
(vs 960 at `-Os`).

### Locked 16 fps, light coverage (frames 193–224)

```
frame                   62.68 ms  37.61 Mcyc  100%
  rs_draw_rings         41.77 ms  25.06 Mcyc   67%
    rs_ring_scan        41.21 ms  24.72 Mcyc   66%  x76     542 us/group
      filter_blend       6.67 ms   4.00 Mcyc   11%  x37807  106 cyc/blend
  rs_buffer_wait        20.86 ms  12.51 Mcyc   33%
  rs_timeline_step      52.4 us   31.5 kcyc     0%
```

Wall min 55.5 / avg 62.7 / max 69.0 ms.

### Trail ramp — the one 2-window window (frames 33–64)

```
frame                  118.76 ms  71.26 Mcyc  100%
  rs_draw_rings         69.03 ms  41.42 Mcyc   58%
    rs_ring_scan        68.35 ms  41.01 Mcyc   58%  x76     899 us/group
      filter_blend      12.80 ms   7.68 Mcyc   11%  x65681  117 cyc/blend
  rs_buffer_wait        49.66 ms  29.80 Mcyc   42%
  rs_timeline_step      65.6 us   39.4 kcyc     0%
```

Wall min 36.6 / avg 118.8 / max 132.4 ms. Peak coverage (65.7 k
blends/frame, 6.3× quadrant) during the trail fill pushes render 6.5 ms
past the window; per-blend scan cost drops to 624 cyc at this density.
Identical in both epochs (the walk is seeded per effect construction).

### Per-pixel figures

Heavy locked window: 46.5 k blends/frame ≈ 4.5× the quadrant at **107 cyc
(0.18 µs) per blend** (162 at `-Os`, 1.51×); scan ≈ 30.9 Mcyc over 46.5 k
blends ⇒ **665 cyc per blended pixel** (960 at `-Os`, 1.44×).

## Column-ISR / DMA marshaling cost

Locked 16 fps windows: `isr_wake` 2.8–3.0%, `isr_pack` 1.35–1.5%,
`isr_dma_submit` 0.21% CPU; the ramp window shows 7.6/6.1/0.21 (2-window
shares plus render-contention inflation). ~1.9 ms ISR per 62.5 ms window in
the locked regime ⇒ ~60.6 ms render budget; the heaviest sustained render
(52 ms) clears it with ~8.5 ms margin — only the ramp peak (69 ms) exceeds
it.

## Summary ranking (heavy locked window, share of the 62.4 ms frame)

1. `rs_ring_scan` (fused `Scan::RingGroup::draw` × 76) — **83%** (51.6 ms)
2. display-window sync (`rs_buffer_wait`) — 16% (10.2 ms)
3. `filter_blend` (inside the scan) — 13% (8.3 ms)
4. everything else — <0.1%

## -O3 vs -Os (both at the fused tip `56d8c854`)

| metric | -Os | -O3 | Δ |
|---|---|---|---|
| `rs_ring_scan`, heavy sustained | 66.8 ms | 51.6 ms | **1.30× faster** |
| scan cyc / blended px | 960 | 665 | 1.44× |
| `filter_blend` per blend | 162 cyc | 107 cyc | 1.51× |
| cadence | ~9.4 fps effective | **~15.6 fps, locked 16** | tier gained |
| image (FLASH / ITCM) | 41.9 / 28.5 KB | 86.2 / 70.4 KB | +106% / +147% |

Vs the pre-fused `-O3` baseline (2026-07-14: in-band scan 54.6 ms,
715 cyc/blend, borderline 16↔8): the fusion bought **−5.6% scan / −7%
per blend at `-O3`** — modest, but exactly the margin that converts the
borderline cadence into a lock. At `-Os` the same change is neutral (see
the twin), so RingSpin is now a prime candidate for the selective-`-O3`
spec (docs/selective_o3_spec.md): its 16 fps tier costs only the fused
scan's translation units at `-O3`.

## Caveats

- Not the shipping config: the full 26-effect roster overflows ITCM at
  `-O3`; single-effect image only.
- All scopes absorb live ISR time; `filter_blend` parents under
  `rs_ring_scan`.
- `-ffast-math` changes float semantics vs the `-Os`/host builds; output
  is not bit-comparable across optimization levels.
- Same fused-path AA-tail divergence note as the `-Os` twin.
- Epoch straddle windows skipped; the ramp window (33–64) includes the
  trail fill.

## Harness

`pio run -e profile_o3 -t upload` with
`PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=RingSpin -D HS_PROFILE_WINDOW=32'`,
then `python tools/profile_capture.py --seconds 150 --out <log>`.
