# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-07-15, **-O3**)

Point-in-time snapshot at the landed fused-scan tip (`7d50b672`); replaces
the 2026-07-14 legacy-path report. Raw capture:
`build/profile_capture_tip_o3_w32.log` (32-frame windows, ~150 s). This is
the global-`-O3` ceiling twin of the shipping
[selective-O3 report](../shipping/profile_displacementfield_teensy_2026-07-15.md);
setup and phase schedule are identical except the optimization level.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env — the -O3 twin of the shipping `-Os` profile: same Phantasm flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE -D HS_PROFILE_WINDOW=32`, but base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` |
| Driver | `POVSegmented<288, 4, 480>` — segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `DisplacementField<288, 144>` at master `7d50b672` (fused stack scan, shipping `thickness = 0.04`) |
| Method | DWT `CYCCNT` scopes + `micros()` wall per `draw_frame`, dumped every 32 frames |
| Reproduce | `pio run -e profile_o3 -t upload` then `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code 89,100 B; ITCM (RAM1 code) 69,816 B; RAM2 free 4,736 B.
(At -Os: 57,236 / 39,416 / 4,736 — FLASH +56%, ITCM +77%.)

**Exactness cross-check** — window frames 385–416: root counter
1,200,381,617 cyc = 2,000,636 µs vs measured wall sum 2,000,639 µs
(Δ ≈ 1.5 ppm).

## Frame cadence

**Every phase runs 1-window cadence: 62.5 ms, 16 fps** — the fused scan
plus `-O3` crosses the dwell and steady-BALLS tiers that the 2026-07-14
legacy `-O3` image (78 ms render) missed. Dwell render is 59.5 ms with
~3 ms of margin; `df_buffer_wait` shrinks to the residual sync idle.

## Phase-by-phase readout (32-frame windows, per-frame averages)

### NOISE dwell (frames 385–416)

```
frame                   62.52 ms  37.51 Mcyc  100%
  df_timeline_step      59.48 ms  35.69 Mcyc   95%
    df_draw_rings       59.38 ms  35.63 Mcyc   95%
      df_fused_scan     45.97 ms  27.58 Mcyc   74%
        filter_blend                                x13615  164 cyc/blend
      df_lut_bake        9.33 ms   5.60 Mcyc   15%  x35.7   261 us/ring
      df_hue_table_prep  1.81 ms   1.09 Mcyc    3%  x22.3    81 us/call
      df_chunk_cull      1.28 ms   0.77 Mcyc    2%  x41.1    31 us/ring
  df_buffer_wait         3.04 ms   1.82 Mcyc    5%
```

Wall min 58.8 / avg 62.5 / max 66.0 ms — **16 fps** (vs 125 ms at -Os).
The second epoch's same window agrees (render 58.3 ms). vs -Os: render
82.6 → 59.5 (**1.39×**), scan 62.7 → 46.0 (1.36×), bake 12.9 → 9.3
(1.38×), hue prep 4.1 → 1.8 (2.2×).

### Steady BALLS (frames 993–1024)

```
frame                   62.46 ms  37.47 Mcyc  100%
  df_timeline_step      47.42 ms  28.45 Mcyc   75%
    df_draw_rings       47.24 ms  28.34 Mcyc   75%
      df_fused_scan     31.52 ms  18.91 Mcyc   50%
        filter_blend                                x12290  165 cyc/blend
      df_lut_bake       10.35 ms   6.21 Mcyc   17%  x41.7   248 us/ring
      df_hue_table_prep  2.46 ms   1.48 Mcyc    4%  x39.8    62 us/call
      df_chunk_cull      1.28 ms   0.77 Mcyc    2%  x41.7    31 us/ring
  df_buffer_wait        15.02 ms   9.01 Mcyc   24%
```

16 fps with 15 ms of idle (vs -Os: 3.8 ms over the window). vs -Os:
render 66.3 → 47.4 (1.40×), scan 43.6 → 31.5 (1.38×).

### Early BALLS (frames 897–928)

```
frame                   62.73 ms  37.64 Mcyc  100%
  df_timeline_step      23.56 ms  14.14 Mcyc   37%
    df_draw_rings       23.45 ms  14.07 Mcyc   37%
      df_fused_scan     19.54 ms  11.73 Mcyc   31%
        filter_blend                                x12208  158 cyc/blend
      df_lut_bake        1.84 ms   1.10 Mcyc    3%  x13.6   135 us/ring
      df_hue_table_prep  0.94 ms   0.56 Mcyc    1%  x11.9
  df_buffer_wait        39.16 ms  23.50 Mcyc   62%
```

### Per-pixel figures

Dwell `filter_blend`: 13,615 blended px/frame at 164 cyc/blend;
`df_fused_scan` ≈ 27.6 Mcyc ⇒ ~2,026 scan cycles per blended pixel
(-Os: ~2,740). Blend itself speeds 194 → 164 cyc (1.18×); the polyline
search gains more (~1.4×) from `-O3` scheduling of the chart arithmetic.

## Column-ISR / DMA marshaling cost

1-window cadence throughout: `isr_wake` ≈ 3.1%, `isr_pack` ≈ 1.5%,
`isr_dma_submit` ≈ 0.2% CPU. All render scopes absorb this.

## Summary ranking (NOISE dwell frame)

1. `df_fused_scan` — 46.0 ms, 74% of frame.
2. `df_lut_bake` — 9.3 ms, 15%.
3. `df_buffer_wait` — 3.0 ms residual sync idle.
4. `df_hue_table_prep` — 1.8 ms, 3%.

**-O3 vs -Os at the fused tip: 1.39× dwell render, and the 8 → 16 fps
cadence crossing for the entire cycle.** This is the measured ceiling for
the selective-`-O3` spec (docs/selective_o3_spec.md): compiling just
`df_fused_scan` + the bake at `-O3` inside the `-Os` shipping image would
recover most of the tier at a fraction of the +30.4 KB ITCM cost.

## Caveats

- Not the shipping config: the full 26-effect roster overflows ITCM at
  `-O3`; single-effect image only.
- All scopes absorb live ISR time; `filter_blend` parents under
  `df_fused_scan`.
- `-ffast-math` changes float semantics vs the `-Os`/host builds; output
  is not bit-comparable across optimization levels.
- Same thickness note as the `-Os` twin: 0.04 here vs 0.035 in all
  2026-07-14 reports.

## Harness

`pio run -e profile_o3 -t upload` with
`PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=DisplacementField -D HS_PROFILE_WINDOW=32'`,
then `python tools/profile_capture.py --seconds 150 --out <log>`.
