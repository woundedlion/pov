# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingSpin`).
Raw capture: `build/prof/ringspin_ship.log`. Third and final RingSpin report of
2026-07-25: the morning sweep caught the old renderer red (64.44 ms, 6/1088
spilled); the row-local RingGroup walk (`2a50eb82`) landed and took it green; a
swept-band trail rework (`83206a52`) briefly took peak to 36 ms but was
**reverted on owner look-rejection** (`99f1dfa4` — the continuous fill
destroyed the discrete-ring trail texture). This report describes the shipping
state at tip `99f1dfa4`: the original look, rendered by the row-local walk.

The roster index ([`../README.md`](../README.md)) ranks this row from a
2026-07-26 11:44 log that postdates this report: the index peak is the current
figure and the numbers below are the earlier capture. `just profile RingSpin`
regenerates this report against a fresh log.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os`, newlib-nano, `USE_DMA_LEDS`, `HS_PROFILE_ENABLE`, `tools/phantasm.ld`, N=4. `RingGroup::draw` now walks its rows inside its own `HS_O3` region (stack interval buffers, row-hoisted slot gate) instead of calling the -Os `scan_region`, whose optimize-attribute mismatch forced the O3 pixel body out of line — one function call per covered pixel in the old image. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, tip `99f1dfa4` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingSpin profile 70 32` |

Image size: `FLASH: code:50796, data:537692, headers:8496` / `RAM1:
variables:314400, code:35976, padding:29560, free:144352` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 929–960 root counter cyc ÷ 600 MHz matches
the measured wall sum within **2.3 ppm**
(`tools/parse_profile.py build/prof/ringspin_ship.log validate`).

## Frame cadence

**Pass aggregate**: peak frame render **56.08 ms** (frames 897–928), spilled
**0/1088 frames (0%)** — 🟢, ~6.4 ms under the 62.5 ms window at the same
overlap crest that pierced it in the morning capture (64.44 ms on the old
renderer, same trajectory — per-window blend counts match). RingSpin still
swings hard within a window (peak window wall min/avg/max ≈ 46/73/126 ms), but
the crest now clears.

## Phase-by-phase readout

Phase schedule: RingSpin is not a cycler; it continuously spins 4 great-circle
ring trails. Single regime; the peak window is below.

### Peak window (frames 897–928)

```
frame                61.59 ms  36.96 Mcyc  100%
  rs_draw_rings      41.43 ms  24.86 Mcyc   67%
    rs_ring_scan     40.80 ms  24.48 Mcyc   98%  x76.0  537 us/group
      filter_blend    8.78 ms   5.27 Mcyc   21%  x52923  100 cyc/blend
  rs_timeline_step    0.09 ms   0.06 Mcyc    0%
  rs_buffer_wait     20.07 ms  12.04 Mcyc   32%
```

Wall min/avg/max = 48.22/61.59/74.80 ms; percentages are of the parent scope.
76 fused ring-group scans per frame (19 trail frames × 4 rings, fused ≤4 per
scan); 52,923 blended px/frame at the crest = **5.10× overdraw** — unchanged
from the old renderer by design (the look depends on it). What changed is the
cost per evaluated pixel: pass-wide `rs_ring_scan` is **32.82 ms/f at 470
cyc/blend**, versus 37.89 ms/f at 541 before the walk landed.

### Per-pixel figures

Quadrant = 10,368 px; `filter_blend` = 52,923 blended px/frame at crest ⇒
**5.10× coverage** (ring overlap; the trail's stacked-ring texture).
`filter_blend` 100 cyc/blend; whole-scan **470 cyc per blended pixel**
pass-wide — below the global-O3 ceiling's 477.

## Column-ISR / DMA marshaling cost

```
isr_wake        1136/frame  min/avg/max 0.64/1.75/11.71 us  cpu 3.23%
isr_pack         142/frame  min/avg/max 6.11/6.61/9.51 us  cpu 1.52%
isr_dma_submit   142/frame  min/avg/max 0.68/0.94/1.06 us  cpu 0.21%
```

(window frames 897–928, the peak window.) Total ISR CPU share **4.96%**.

## Summary ranking

1. `rs_ring_scan` — 76 fused group scans over 54k blended pixels at crest; the
   entire render cost. Now fully in-region: no per-pixel call, no arena/CSG row
   machinery, band gate hoisted per row.
2. `rs_buffer_wait` — display-sync idle.

**The shipping image now beats the global-O3 ceiling on this path** — 470 vs
477 cyc/blend, peaks 56.08 vs 54.45 ms (1.03×) — because the bespoke walk
removes structural overhead the ceiling's general `scan_region` still pays
(arena scratch + memset, interval sort/normalize calls, udiv buffer indexing).
The remaining crest cost is the 5.25× overdraw itself, which is the effect's
visual identity: a swept-band rework that collapsed it to 36 ms was reverted
because the continuous fill destroyed the discrete-ring texture (see
`99f1dfa4`). Further optimization must preserve per-frame ring striation.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `rs_ring_scan`; its calls ≈ blended pixels.
- Output is coverage-identical to the pre-walk renderer (same pixels, same
  blend order); only the interval plumbing changed.
- Crest peaks remain capture-dependent (which rotation phase the 70 s samples);
  with the walk the crest tops out ~6 ms under the window instead of ~2 ms
  over, so the colour is no longer a coin flip.
- Working tree tip `99f1dfa4`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingSpin`,
`HS_PROFILE_WINDOW=32`; `just profile RingSpin` = build + flash + capture.
