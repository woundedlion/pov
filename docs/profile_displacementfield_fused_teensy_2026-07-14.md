# DisplacementField (fused scan) on-device profile — Teensy 4.0, segmented mode (2026-07-14)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture:
`build/profile_capture_fused_w32.log` (32-frame windows, ~150 s). Delta
capture against the same-day legacy baseline
[profile_displacementfield_teensy_2026-07-14.md](profiles/Os/profile_displacementfield_teensy_2026-07-14.md);
Setup/ISR context identical unless noted.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`) + `-D HS_PROFILE_ENABLE -D HS_PROFILE_WINDOW=32` |
| Driver | `POVSegmented<288, 4, 480>`, this board = segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `DisplacementField<288, 144>`, **fused single-pass scan** (`Scan::DistortedRingStack`) |
| Method | DWT `CYCCNT` scopes (`HS_PROFILE`) + `micros()` wall per `draw_frame`, dumped every 32 frames |
| Reproduce | `just profile DisplacementField` |

**Exactness cross-check** — window frames 385–416: root counter
2,401,138,148 cyc = 4,001,897 µs vs measured wall sum 4,001,900 µs
(Δ ≈ 0.75 ppm).

## Frame cadence

Same 62.5 ms display-window quantization as the baseline (one quadrant ≈
10,370 px per `draw_frame`). Observed:

- **NOISE dwell / steady BALLS**: still 2-window, **125 ms (8 fps)** — but
  steady BALLS render is now 64.4 ms, only **1.9 ms above the window**.
- **Fade-in ramp, fade-out, transition, early BALLS**: **62.5 ms (16 fps)**.
  The ramp windows (65–96) ran dwell-priced (125 ms) on the legacy path —
  a full cadence tier recovered there.

## Phase-by-phase readout (32-frame windows, per-frame averages)

Same phase schedule as the baseline (fade-in ≈ 0–150, dwell ≈ 150–750,
fade-out ≈ 750–900, BALLS ≈ 900+).

### NOISE dwell (frames 385–416) — vs legacy in brackets

```
frame                  125.06 ms  75.04 Mcyc  100%
  df_timeline_step      79.81 ms  47.88 Mcyc   64%   [105.51  84%]
    df_draw_rings       79.69 ms  47.81 Mcyc   64%   [105.41]
      df_fused_scan     59.14 ms  35.48 Mcyc   47%   [df_ring_scan 86.79  69%]
      df_lut_bake       13.69 ms   8.21 Mcyc   11%  x35.7   383 us/ring
      df_hue_table_prep  4.06 ms   2.44 Mcyc    3%  x22.1   184 us/call
      df_chunk_cull      1.48 ms   0.89 Mcyc    1%  x41.9    35 us/ring
  df_buffer_wait        45.25 ms  27.15 Mcyc   36%   [19.56  16%]
```

Wall min 124.6 / avg 125.1 / max 126.0 ms. The fused scan is **−31.9%**
(86.79 → 59.14 ms) and total render **−24.4%** (105.5 → 79.8 ms); ring
counts and bake cost match the baseline (same rings survive the same
culls — the bake is untouched). The saving is the per-ring recompute of
the shared-axis pixel frame (dot/fast_acos/fast_atan2 across the ~7×
overlapping ring bands), now paid once per pixel.

### Steady BALLS (frames 993–1024)

```
frame                  125.05 ms  75.03 Mcyc  100%
  df_timeline_step      64.37 ms  38.62 Mcyc   51%   [79.50  64%]
    df_draw_rings       64.13 ms  38.48 Mcyc   51%
      df_fused_scan     40.39 ms  24.24 Mcyc   32%   [df_ring_scan 58.42  47%]
      df_lut_bake       14.70 ms   8.82 Mcyc   12%  x44.8   328 us/ring
      df_hue_table_prep  5.34 ms   3.20 Mcyc    4%  x40.1   133 us/call
      df_chunk_cull      1.55 ms   0.93 Mcyc    1%  x44.9    35 us/ring
  df_buffer_wait        60.66 ms  36.40 Mcyc   48%
```

Scan −30.8%, render 79.5 → 64.4 ms — **1.9 ms short of the 62.5 ms
window**. The bake (14.7 ms) is now the swing item: the ledger's
candidate #4 (per-(ring,ball) azimuth-interval bake; the noise field is
exactly zero in this phase) would buy the tier.

### Early BALLS (frames 897–928)

```
frame                   62.89 ms  37.73 Mcyc  100%
  df_timeline_step      36.62 ms  21.97 Mcyc   58%   [34.68  55%]
    df_draw_rings       36.48 ms  21.89 Mcyc   58%
      df_fused_scan     30.34 ms  18.20 Mcyc   48%
      df_lut_bake        2.82 ms   1.69 Mcyc    4%  x13.5
      df_hue_table_prep  1.99 ms   1.19 Mcyc    3%  x12.0
  df_buffer_wait        26.27 ms  15.76 Mcyc   42%
```

16 fps, unchanged tier, but render is **+5.5%** vs legacy (34.7 → 36.6):
the ~25 flat rings/frame that took the `-Os` explicit-SDF `draw_flat`
path now ride the shared fused pass as zero-knot candidates. Cadence is
unaffected; a flat fast path inside the candidate loop could claw this
back if it ever matters.

### Ramp-in (frames 65–96)

```
frame                   62.79 ms  37.67 Mcyc  100%
  df_timeline_step      57.17 ms  34.30 Mcyc   91%   [legacy: dwell-priced]
    df_fused_scan       42.26 ms  25.36 Mcyc   67%
    df_lut_bake         11.24 ms   6.74 Mcyc   18%
  df_buffer_wait         5.62 ms   3.37 Mcyc    9%
```

**16 fps where the legacy path was already at 125 ms** — the fused scan
holds the ramp under one window until the amplitude nears full.

### Per-pixel figures

Dwell `filter_blend`: 12,020 blends/frame at 191 cyc/blend (≈2.3 Mcyc);
`df_fused_scan` ≈ 35.5 Mcyc over ~10.4 k quadrant pixels ⇒ ~3,400 cyc per
covered pixel (was ~4,780) — the remainder is the per-candidate polyline
search, now the dominant irreducible term.

## Column-ISR / DMA marshaling cost

Unchanged from the baseline (same driver, same windows): `isr_wake` ≈
4.7%, `isr_pack` ≈ 3.0%, `isr_dma_submit` ≈ 0.2% CPU in 2-window dwell;
halved shares in 1-window phases. All render scopes absorb this.

## Summary ranking (NOISE dwell frame)

1. `df_fused_scan` — 59.1 ms, 47% of frame (**was 86.8 ms / 69%**).
2. `df_buffer_wait` — 45.3 ms display sync idle (grows as render shrinks).
3. `df_lut_bake` — 13.7 ms, 11%.
4. `df_hue_table_prep` — 4.1 ms, 3%.

Host A/B (native -O2, interleaved, FNV bit-identical over 2000 frames incl.
sabotage check): dwell 8.50 → 7.04 ms/frame (−17%), balls 5.08 → 4.28
(−16.5%). Device gains are larger (−31%) — the in-order M7 pays more for
the redundant per-ring trig than the desktop OoO core.

## Caveats

- All scopes absorb live ISR time; `filter_blend` parents under
  `df_fused_scan` and counts all blended pixels.
- Captured on the tree that landed the fused scan
  (`Scan::DistortedRingStack` + `SDF::DistortedRing::distance_from_frame`)
  and the `suppress_pole_fill` pole-row fix.
- `-Os` shipping config; epoch (120 s) reconstructs the effect ~frame 1030,
  so later windows straddle.

## Harness

`just profile DisplacementField [seconds]`; knobs via
`PLATFORMIO_BUILD_FLAGS` (`HS_PROFILE_TARGET`, `HS_PROFILE_WINDOW`).
