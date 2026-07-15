# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-07-15, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `build/profile_displacementfield_selo3_v2.log`; same-tip
baselines: `build/profile_displacementfield_base_os.log` (pre-region `-Os`),
`build/profile_displacementfield_base_o3.log` (global-O3 twin). The `Os/`
and `O3/` 2026-07-15 reports remain the per-level references for the
pre-region tree.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: shipping Phantasm flags (`-Os`, newlib-nano, DMA LEDs, N=4) + `HS_PROFILE_ENABLE` + the landed `HS_O3` regions (blend sink, `Scan::DistortedRingStack`, sdf.h `HS_O3_FN`s, `draw_rings`) |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture (full NOISE→BALLS cycle + epoch wrap) |
| Reproduce | `just profile DisplacementField` (the `-Os` `profile` env measures the shipping selective-O3 config) |

Image size: `FLASH: code:62492, data:79068, headers:8960` / `RAM1:
variables:345824, code:42904, padding:22632, free:112928` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 193–224 root counter
2,399,804,966 cyc ÷ 600 MHz = 3,999,675 µs vs measured wall sum
3,999,682 µs — **1.8 ppm**.

## Frame cadence

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Observed cadence across the cycle: **BALLS phase, the NOISE ramps, and the
steady NOISE regime lock 16 fps** (62.1–62.9 ms window averages — all of
this was 8 fps at `-Os` through the dwell); the **peak NOISE dwell**
(~14 of 55 windows) still averages 125 ms (8 fps) — its render is
63–66 ms, just over one window. The global-O3 ceiling is itself marginal in
that stretch (render 59.7 ms, window average 62.56 ms, single frames
spilling to 65 ms). `df_buffer_wait` is the display-sync round-up idle.

## Phase-by-phase readout

Phase schedule: NOISE fade-in → NOISE dwell (peak) → fade-out → BALLS →
repeat; the palette wipes every ~3 s throughout.

### NOISE dwell, peak (window frames 193–224)

```
frame                  124.99 ms  75.00 Mcyc  100%
  df_timeline_step      66.41 ms  39.85 Mcyc   53%
    df_draw_rings       66.28 ms  39.77 Mcyc   53%
      df_hue_table_prep  3.93 ms   2.36 Mcyc    3%  x22.6  174 us/call
      df_lut_bake       11.73 ms   7.04 Mcyc    9%  x39.6  296 us/call
      df_chunk_cull      1.31 ms   0.78 Mcyc    1%  x44.3   29 us/call
      df_fused_scan     48.16 ms  28.89 Mcyc   39%  x1     48.2 ms/scan
        filter_blend     3.75 ms   2.25 Mcyc    3%  x13642 165 cyc/blend
  df_buffer_wait        58.58 ms  35.15 Mcyc   47%
  df_prepare_fields      0.13 us     84 cyc     0%
```

Wall min/avg/max = 124.2/125.0/125.7 ms. Render (draw_rings) is 66.3 ms —
3.8 ms over the window, so the frame takes 2 windows and the idle balloons.
The same window at `-Os` was scan 60–62 + bake 16.0 + hue 4.6 (render
~78 ms); at global -O3, scan 44.9 + bake 10.6 + hue 1.7 (render 59.7 ms,
marginal lock). The residual 3–6 ms vs the ceiling is the still-`-Os` OKLab
hue chain and the `FastNoiseLite`-driven bake interior (measured; see the
spec ledger — a hue-member-only promotion was measured dead and reverted).

### NOISE steady (window frames 1089–1120)

```
frame                   62.51 ms  37.51 Mcyc  100%
  df_timeline_step      46.76 ms  28.06 Mcyc   75%
    df_draw_rings       46.52 ms  27.91 Mcyc   74%
      df_hue_table_prep  2.55 ms   1.53 Mcyc    4%  x33.5   76 us/call
      df_lut_bake       10.40 ms   6.24 Mcyc   17%  x39.8  261 us/call
      df_chunk_cull      1.25 ms   0.75 Mcyc    2%  x43.8   29 us/call
      df_fused_scan     30.72 ms  18.43 Mcyc   49%  x1     30.7 ms/scan
        filter_blend     3.22 ms   1.93 Mcyc    5%  x12295 157 cyc/blend
  df_buffer_wait        15.73 ms   9.44 Mcyc   25%
  df_prepare_fields      9.2 us    5.5 kcyc     0%
```

Wall min/avg/max = 58.9/62.5/66.0 ms — **16 fps locked** with 25% idle. At
`-Os` this regime ran 2 windows; the selective regions carry the crossing.

### BALLS phase (window frames 897–928)

```
frame                   62.88 ms  37.73 Mcyc  100%
  df_timeline_step      24.16 ms  14.50 Mcyc   38%
    df_draw_rings       24.02 ms  14.41 Mcyc   38%
      df_hue_table_prep  1.99 ms   1.19 Mcyc    3%  x12.0  165 us/call
      df_lut_bake        1.77 ms   1.06 Mcyc    3%  x13.6  131 us/call
      df_chunk_cull      0.38 ms   0.23 Mcyc    1%  x14.0   27 us/call
      df_fused_scan     19.10 ms  11.46 Mcyc   31%  x1     19.1 ms/scan
        filter_blend     3.06 ms   1.84 Mcyc    5%  x12208 150 cyc/blend
  df_buffer_wait        38.72 ms  23.23 Mcyc   61%
  df_prepare_fields      1.4 us    859 cyc      0%
```

Wall min/avg/max = 52.5/62.9/68.9 ms — 16 fps with wide idle; the ball
phase was already light and rides along.

### Per-pixel figures

Dwell peak: 13,642 blended px/frame vs the 10,368 px quadrant = 1.3×
coverage; 165 cyc/blend (DF's shader does per-fragment hue-table sampling —
heavier than RingSpin's 84); fused-scan cost ≈ 2,118 cyc per blended pixel
(the stack's candidate loop dominates, not the blend).

## Column-ISR / DMA marshaling cost

```
isr_wake        2303/frame  min/avg/max 0.64/2.9/151 us   cpu 5.36%
isr_pack         288/frame  min/avg/max 4.3/15.8/149 us   cpu 3.63%
isr_dma_submit   288/frame  min/avg/max 0.61/0.97/7.5 us  cpu 0.22%
```

(2-window dwell frames double the per-frame ISR counts.) ISR CPU share ≈
9.2% in the dwell ⇒ ~56.7 ms of render budget per 62.5 ms window — the peak
dwell's 66.3 ms render misses it; the steady regime's 46.5 ms fits.

## Summary ranking

1. `df_fused_scan` — 31–39% of the (multi-window) frame, 19.1–48.4 ms:
   within ~4% of the global-O3 ceiling in the dwell (45.7 vs 44.0 ms
   medians).
2. `df_lut_bake` — 9–17%, 10.4–11.7 ms in NOISE: ~1.2 ms over the ceiling
   figure (FastNoiseLite interior still `-Os`).
3. `df_hue_table_prep` — 3–4%: ~2 ms over the ceiling (OKLab chain still
   `-Os`).
4. `filter_blend` — 3–5%, 150–165 cyc/blend at the -O3 sink (190 at `-Os`).

Selective -O3 closes ~88% of the dwell's -Os → -O3 hot-scope gap and
crosses everything except the peak-dwell stretch to 16 fps.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `df_fused_scan`; calls ≈ blended pixels.
- Epoch-boundary windows (frame calls > 32) were skipped for the
  representative picks.
- Shipping config is `-Os` + selective `HS_O3` regions — this is the
  shipping image's profile.
- Captured at the DF draw-path region tip (`e98f4651`); the later hue-member
  promotion was confirmed no-change (`build/..._selo3_v3.log`) and reverted
  (`7284bb09`), so this capture matches the final tree.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`,
`HS_PROFILE_WINDOW=32`; `just profile DisplacementField` = build + flash +
capture.
