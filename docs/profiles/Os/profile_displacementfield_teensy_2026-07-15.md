# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-07-15)

Point-in-time snapshot at the landed fused-scan tip (`7d50b672`); replaces
the 2026-07-14 legacy-path report. Regenerate with
`just profile DisplacementField`. Raw captures:
`build/profile_capture_tip_w32.log` (32-frame windows, ~150 s),
`build/profile_capture_tip_w128.log` (128-frame windows, ~260 s, two
epochs), plus the A/B decomposition captures (`..._preopt2_w32.log`,
`..._hoist_w32.log`, `..._opt2_w32.log`). `-O3` twin:
[../O3/profile_displacementfield_teensy_2026-07-15.md](../O3/profile_displacementfield_teensy_2026-07-15.md).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile` env = Phantasm shipping flags (`-Os`, newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE` (`HS_PROFILE_WINDOW=32` for the fine pass) |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as segment 0 (master), flywheel + DMA LED ISRs live |
| Effect | `DisplacementField<288, 144>` at master `7d50b672`: fused single-pass scan (`Scan::DistortedRingStack`), flat fast path, chart-sin hoist, shipping `thickness = 0.04` |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped per window then reset |
| Reproduce | `just profile DisplacementField` |

Image size: FLASH code 57,236 B; ITCM (RAM1 code) 39,416 B; RAM2 free 4,736 B.

**Exactness cross-check** — window frames 385–416 (w32): root counter
2,401,109,540 cyc = 4,001,849 µs vs measured `micros()` window sum
4,001,860 µs (Δ ≈ 2.7 ppm).

**Thickness note:** the 2026-07-14 captures ran the profiling tree's 0.035
default; this tree runs the landed `params.thickness = 0.04f`, which alone
costs ~+10% scan. Cross-date comparisons are invalid; the same-thickness
A/Bs are in the Summary section.

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; one
`draw_frame` renders one quadrant (72-row band × 144-column half ≈
10,370 px). Wall time quantizes to whole windows:

- **NOISE dwell / steady BALLS**: 2-window cadence, **125 ms (8 fps)**.
- **Fade-in, fade-out, early BALLS**: **62.5 ms (16 fps)**, mixing 2-window
  frames near the phase edges (fade-out min 57.5 / max 125.3 ms).

`df_buffer_wait` (the `buffer_free()` spin inside the Canvas ctor) is the
round-up idle to the next window boundary; render numbers below exclude it.

## Phase-by-phase readout (32-frame windows, per-frame averages)

Phase schedule: fade-in ≈ frames 0–150, NOISE dwell ≈ 150–750, fade-out ≈
750–900, BALLS from ≈ 900; epoch (120 s) reconstructs the effect ≈ frame
1030.

### NOISE dwell — the hot phase (frames 385–416)

```
frame                  125.06 ms  75.03 Mcyc  100%
  df_timeline_step      82.60 ms  49.56 Mcyc   66%
    df_draw_rings       82.48 ms  49.49 Mcyc   66%
      df_fused_scan     62.70 ms  37.62 Mcyc   50%
        filter_blend                                x13748  194 cyc/blend
      df_lut_bake       12.90 ms   7.74 Mcyc   10%  x35.8   360 us/ring
      df_hue_table_prep  4.06 ms   2.44 Mcyc    3%  x22.1   184 us/call
      df_chunk_cull      1.50 ms   0.90 Mcyc    1%  x42.0    36 us/ring
  df_buffer_wait        42.45 ms  25.47 Mcyc   34%
  df_prepare_fields      0.16 us      98 cyc    0%
```

Wall min 124.3 / avg 125.1 / max 125.9 ms. ~35.8 of 48 rings bake + scan.
The w128 pass puts full-dwell windows (257–384) at 79.2 (epoch 1) and
77.8 ms (epoch 2) render — dwell cost breathes ~±5% with noise content but
the two epochs agree.

### Steady BALLS (frames 993–1024)

```
frame                  125.06 ms  75.04 Mcyc  100%
  df_timeline_step      66.28 ms  39.77 Mcyc   53%
    df_draw_rings       66.04 ms  39.63 Mcyc   53%
      df_fused_scan     43.59 ms  26.16 Mcyc   35%
        filter_blend                                x12301  196 cyc/blend
      df_lut_bake       13.48 ms   8.09 Mcyc   11%  x44.8   301 us/ring
      df_hue_table_prep  5.27 ms   3.16 Mcyc    4%  x40.1   132 us/call
      df_chunk_cull      1.58 ms   0.95 Mcyc    1%  x44.9    35 us/ring
  df_buffer_wait        58.76 ms  35.26 Mcyc   47%
```

Render 66.3 ms — **3.8 ms above the 62.5 ms window**. The bake (13.5 ms)
is the swing item: the ledger's candidate #4 (per-(ring,ball)
azimuth-interval bake; the noise field is exactly zero in this phase)
would buy the 16 fps tier.

### Early BALLS (frames 897–928)

```
frame                   62.95 ms  37.77 Mcyc  100%
  df_timeline_step      33.88 ms  20.33 Mcyc   54%
    df_draw_rings       33.74 ms  20.24 Mcyc   54%
      df_fused_scan     28.17 ms  16.90 Mcyc   45%
        filter_blend                                x12203  182 cyc/blend
      df_lut_bake        2.32 ms   1.39 Mcyc    4%  x13.4   173 us/ring
      df_hue_table_prep  1.94 ms   1.16 Mcyc    3%  x11.9
      df_chunk_cull      0.44 ms   0.26 Mcyc    1%  x14.3
  df_buffer_wait        29.07 ms  17.44 Mcyc   46%
```

16 fps (wall min 48.4 / avg 62.95 / max 71.7). The flat fast path keeps
the many undisplaced rings at par with the old explicit flat SDF: render
matches the pre-fused shape of this window at the same thickness.

### Fades

- **Fade-in (frames 1–32)**: global minimum, wall min 19.8 ms; render
  climbs with the amplitude ramp (window avg 42.8 ms render, 16 fps).
- **Fade-out (frames 833–864)**: mixed cadence (57.5–125.3 ms), render
  62.4 ms average falling as the band narrows.

### Per-pixel figures

Dwell `filter_blend`: 13,748 blended px/frame (the 48-ring stack at 0.04
thickness over-covers the 10,368-px quadrant) at 194 cyc/blend;
`df_fused_scan` ≈ 37.6 Mcyc ⇒ ~2,740 scan cycles per blended pixel
(0.035-era fused: ~3,400; legacy per-ring: ~4,780). The residue is the
per-candidate polyline search — the dominant irreducible term.

## Column-ISR / DMA marshaling cost

Dwell windows (w32, per-frame rates): `isr_wake` 2,304/frame, avg
1,551 cyc, **4.8% CPU**; `isr_pack` 288/frame, avg 7,964 cyc, **3.1%
CPU**; `isr_dma_submit` 288/frame, avg 580 cyc, **0.2% CPU**. One-window
phases halve the shares (3.1 / 1.4 / 0.2%).

- ISR overhead ≈ 8% CPU in dwell ⇒ ~57.5 ms of render budget per 62.5 ms
  window.
- Dwell render 82.6 ms needs **−30%** for 16 fps; steady BALLS 66.3 ms
  needs **−6%**.

## Summary ranking (NOISE dwell frame)

1. `df_fused_scan` — 62.7 ms, 50% of frame.
2. `df_buffer_wait` — 42.5 ms display-sync idle.
3. `df_lut_bake` — 12.9 ms, 10%.
4. `df_hue_table_prep` — 4.1 ms, 3%.

Same-thickness A/B decomposition (dwell / steady-BALLS / early-BALLS
scan, ms/frame):

| image | scan |
|---|---|
| d463fe63 fused + flat fast path | 64.9 / 45.5 / 28.3 |
| 98dfcdc5 + stack-level gap test + sin hoist | 66.9 / 48.2 / 33.9 |
| **7d50b672 sin hoist only (landed)** | **62.7 / 43.6 / 28.2** |

The stack-level 3-chunk gap test was a net device tax (early-BALLS +20%)
despite host −4% — post-reject work is thin and `-Os` on the in-order
core pays for per-candidate branching; it was reverted. The chart-sin
hoist (−3.3% / −4.2%, bit-identical) landed. Remaining levers: the
selective `-O3` spec (docs/selective_o3_spec.md) for the dwell's 30% —
the `-O3` twin already runs the whole cycle at 16 fps — and ledger
candidate #4 for the steady-BALLS 6%.

## Caveats

- All scopes absorb live ISR time (CYCCNT free-runs); `filter_blend`
  parents under `df_fused_scan` and counts every blended pixel.
- Epoch (120 s) reconstructs the effect ≈ frame 1030; later windows
  straddle the teardown tail.
- `-Os` shipping config; see the `-O3` twin for the optimization-level
  cost.
- A parallel session flashed a different profile image mid-run once; the
  committed captures were verified `effect=DisplacementField` in every
  window header.

## Harness

`just profile DisplacementField [seconds]`; knobs via
`PLATFORMIO_BUILD_FLAGS` (`HS_PROFILE_TARGET`, `HS_PROFILE_WINDOW`).
