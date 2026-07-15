# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-07-15)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw captures: `build/profile_capture_tip_w32.log` (landed tip, 32-frame
windows, ~150 s), plus the A/B decomposition captures
(`..._preopt2_w32.log`, `..._hoist_w32.log`, `..._opt2_w32.log`). Follows
[profile_displacementfield_fused_teensy_2026-07-14.md](profile_displacementfield_fused_teensy_2026-07-14.md);
Setup/ISR context identical. **Thickness note:** all 2026-07-14 captures ran
the profiling tree's 0.035 default; these run the shipping
`params.thickness = 0.04f`, which alone costs ~+10% scan — same-day deltas
below are same-thickness A/Bs.

## Frame cadence

- **NOISE dwell / steady BALLS**: 2-window, 125 ms (8 fps). Steady BALLS
  render is 66.3 ms — 3.8 ms above the 62.5 ms window.
- **Early BALLS / fade edges**: 62.5 ms (16 fps); ramp-in mixes cadences as
  the amplitude climbs (0.04 crosses the window earlier than 0.035 did).

## Phase readout (32-frame windows, per-frame averages, landed tip 7d50b672)

### NOISE dwell (frames 385–416)

```
frame                  125.06 ms  75.04 Mcyc  100%
  df_timeline_step      82.60 ms  49.56 Mcyc   66%
    df_fused_scan       62.70 ms  37.62 Mcyc   50%
    df_lut_bake         12.90 ms   7.74 Mcyc   10%  x35.7
    df_hue_table_prep    4.06 ms   2.44 Mcyc    3%  x22.1
  df_buffer_wait        42.45 ms  25.47 Mcyc   34%
```

### Steady BALLS (frames 993–1024)

```
frame                  125.06 ms  75.06 Mcyc  100%
  df_timeline_step      66.28 ms  39.77 Mcyc   53%
    df_fused_scan       43.59 ms  26.16 Mcyc   35%
    df_lut_bake         13.48 ms   8.09 Mcyc   11%  x44.8
    df_hue_table_prep    5.27 ms   3.16 Mcyc    4%  x40.1
  df_buffer_wait        58.76 ms  35.26 Mcyc   47%
```

### Early BALLS (frames 897–928)

```
frame                   62.95 ms  37.77 Mcyc  100%
  df_timeline_step      33.88 ms  20.33 Mcyc   54%
    df_fused_scan       28.17 ms  16.90 Mcyc   45%
  df_buffer_wait        29.07 ms  17.44 Mcyc   46%
```

## Same-thickness A/B decomposition (dwell / steady-BALLS / early-BALLS scan)

| image | scan ms/frame |
|---|---|
| d463fe63 (fused + flat fast path) | 64.9 / 45.5 / 28.3 |
| 98dfcdc5 (+ stack-level gap test + sin hoist) | 66.9 / 48.2 / 33.9 |
| 7d50b672 (sin hoist only — **landed**) | **62.7 / 43.6 / 28.2** |

The stack-level 3-chunk gap test was a net device tax (early-BALLS +20%)
despite host −4%, echoing the cos-space prefilter result: post-reject work
is thin, and `-Os` on the in-order core pays for per-candidate branching.
The per-pixel chart-sin hoist alone is −3.3% / −4.2% and bit-identical; it
landed (7d50b672), the gap test did not.

## Standing at the shipping config

Dwell render 82.6 ms vs the 62.5 ms window (target < 55): the fused scan is
search-dominated; remaining levers are the selective `-O3` spec
(docs/selective_o3_spec.md) for scan + bake, and the per-(ring,ball)
azimuth-interval BALLS bake (ledger candidate #4) for the 3.8 ms
steady-BALLS gap.

## Harness

`just profile DisplacementField [seconds]`; knobs via
`PLATFORMIO_BUILD_FLAGS` (`HS_PROFILE_TARGET`, `HS_PROFILE_WINDOW`).
