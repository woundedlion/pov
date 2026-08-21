# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-08-21, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile SphericalHarmonics`).
Raw capture: `build/prof/sphericalharmonics_ship.log`. Replaces the 2026-07-25
report, whose numbers predate both a month of unrelated engine work and the
`Scan::Shader` rewrite measured below.

Captured on **COM3** (two boards attached; `HS_TEENSY_PORT` pinned).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. `Scan::Shader::draw_typed` is `HS_O3_FN`, so the whole per-pixel walk is an O3 region; the shader body inlines into it. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, tip `419426cc4` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 220 s capture, epoch stretched to 2048 revs, ordered 24-mode cycle (`HS_PROFILE_ORDERED_CYCLE`). Modes morph back-to-back, so per-mode rows use the peak window nearest each anchor. |
| Reproduce | `bash tools/profile_one.sh SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` |

Image size: `FLASH: code:38872, data:145976, headers:8688` / `RAM1:
variables:315008, code:20664, padding:12104, free:176512` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 2449–2464 root counter cyc ÷ 600 MHz
matches the measured wall sum within **0.3 ppm**
(`tools/parse_profile.py … validate`, all 8 checks pass).

## Frame cadence

**Pass aggregate** (`buckets`): all 24 modes green — worst peak frame render
**12.29 ms** (mode 21, frames 913–928), spilled **0/3488 frames (0%)** — 🟢.
Peaks span 7.63–12.29 ms.

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Render is ~8–12 ms against that 62.5 ms window, so `canvas_buffer_wait` is
84% of the frame — the round-up idle to the next display flip, by design.
SphericalHarmonics remains one of the cheapest cyclers, with ~5× headroom.

## Phase-by-phase readout

Phase schedule: SphericalHarmonics morphs continuously through 24 modes
(`Mode:` marker), no held plateau; the regime below is the peak window nearest
the mode-21 anchor.

### Mode 21 peak window (frames 913–928)

```
frame            62.45 ms  37.47 Mcyc  100%
  sh_rasterize   12.12 ms   7.27 Mcyc   19%
  sh_timeline_step 0.02 ms  0.01 Mcyc    0%
  canvas_buffer_wait 50.29 ms 30.17 Mcyc  81%
```

Wall min/avg/max = 62.06/62.45/62.72 ms; percentages are of the parent scope.
`sh_rasterize` — now a `Scan::Shader::draw` walk rather than an SDF
rasterization — is the whole render. **`filter_blend` no longer appears
anywhere in the tree.** It was a child of `sh_rasterize` in every prior report
(10,658 calls/frame, 0.53 ms/frame) because the old path plotted each shaded
pixel through the terminal `Pipeline`; `Scan::Shader` assigns to the canvas
instead, so the blend is gone rather than merely cheaper.

### Per-preset table

Cycle wrap confirmed (`cycle wraps to 0`, 24 distinct modes, all visited).
Peak/spilled from `buckets`; peak-window figures per the back-to-back-morph
method note. All 24 modes green, 8 clean windows per mode.

| Mode | Peak render ms | `sh_rasterize` ms/f | fps |
|--:|--:|--:|--:|
| 21 | 12.29 | 12.12 | 16 |
| 20 | 12.29 | 12.11 | 16 |
| 19 | 11.12 | 10.92 | 16 |
| 22 | 11.10 | 10.91 | 16 |
| 13 | 10.62 | 10.44 | 16 |
| 12 | 10.61 | 10.43 | 16 |
| 18 | 9.68 | 9.53 | 16 |
| 23 | 9.66 | 9.52 | 16 |
| 14 | 9.31 | 9.12 | 16 |
| 11 | 9.31 | 9.11 | 16 |
| 24 | 8.93 | 8.77 | 16 |
| 17 | 8.93 | 8.77 | 16 |
| 7 | 8.81 | 8.61 | 16 |
| 6 | 8.81 | 8.61 | 16 |
| 16 | 8.63 | 8.44 | 16 |
| 10 | 8.52 | 8.37 | 16 |
| 15 | 8.49 | 8.36 | 16 |
| 9 | 8.22 | 8.06 | 16 |
| 1 | 8.20 | 8.02 | 16 |
| 8 | 8.12 | 7.96 | 16 |
| 5 | 8.08 | 7.95 | 16 |
| 4 | 7.82 | 7.64 | 16 |
| 3 | 7.66 | 7.47 | 16 |
| 2 | 7.63 | 7.46 | 16 |

Cost tracks harmonic degree: the top modes are the l=4 band, the cheapest the
low-l ones, because `SHMath::spherical_harmonic` is the per-pixel work.

### A/B against the pre-rewrite path

Matched pair, captured back-to-back on the **same board** with the **same tip**
(`419426cc4`) differing only in `effects/SphericalHarmonics.h`, reverted to its
`4fbf81e70^` content for the "before" side from a scratch worktree
(`HS_PROFILE_TREE`). So nothing but the shading path varies.

| | before (SDF path) | after (`Scan::Shader`) | Δ |
|---|--:|--:|--:|
| Worst mode peak render | 14.49 ms | 12.29 ms | **−15.2%** |
| Worst mode `sh_rasterize` | 14.30 ms/f | 12.12 ms/f | −15.2% |
| Cheapest mode peak render | 9.86 ms | 7.63 ms | **−22.6%** |
| Peak span across 24 modes | 9.86–14.49 ms | 7.63–12.29 ms | — |
| `filter_blend` | 10,658 calls/f, 0.53 ms/f | absent | removed |
| Spilled | 0/3488 | 0/3488 | — |
| `profile` FLASH code | 39,736 B | 38,872 B | −864 B |
| `profile` RAM1 code (ITCM) | 21,528 B | 20,664 B | −864 B |

The absolute saving is **~2.2 ms/frame on every mode**, which is the expected
shape: the removed work (`DistanceResult` round-trip, `Fragment` scratch,
coverage/threshold branches, type-erased fragment-shader call, terminal
`Pipeline` blend) is a constant per-pixel overhead over a constant pixel count,
while the harmonic evaluation that varies by mode is untouched. The percentage
is therefore largest on the cheapest modes.

Renders are bit-identical across the change (verified natively: FNV-1a folds
over every displayed frame at 288×144 and 96×20 over 200 frames, with Debug BB
on and at amplitude 0.2 and 7.0), so this is pure overhead removal.

### Per-pixel figures

Quadrant = 10,368 px, and the shader writes every pixel of the clip band.
Worst mode: 7.27 Mcyc ÷ 10,368 px = **701 cyc per pixel**, against 828 cyc per
pixel before the rewrite — **127 cyc/px** of per-pixel scaffolding removed.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.70/1.95/11.29 us  cpu 3.59%
isr_pack         144/frame  min/avg/max 6.23/6.82/9.02 us  cpu 1.57%
isr_dma_submit   144/frame  min/avg/max 0.68/0.93/1.06 us  cpu 0.21%
```

Total ISR CPU share **5.37%**, unchanged by the rewrite (it is driver-side).

## Summary ranking

1. `sh_rasterize` — 19% of the frame, 12.12 ms at the worst mode: the
   full-sphere harmonic shade. Now the *only* render scope.
2. `sh_timeline_step` — 0.02 ms, negligible.

`filter_blend` has left the ranking entirely; it was second at 1.5 ms in the
2026-07-25 report.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- There is no `filter_blend` row to parent: `Scan::Shader` assigns to the
  canvas and takes no pipeline, so no plot-time filter stage runs. The effect
  carries no filters, so nothing is lost — an effect needing the filter chain
  could not use this path.
- Back-to-back morphs: per-mode rows are peak windows nearest each anchor.
- Dwell knobs (`HS_PROFILE_ORDERED_CYCLE`, `HS_PROFILE_EPOCH_REVS=2048`) change
  how long a mode holds, not per-frame cost.
- The "before" side of the A/B ran with an uncommitted revert of
  `effects/SphericalHarmonics.h` in a scratch worktree; the "after" side is
  clean at `419426cc4`.
- SphericalHarmonics no longer participates in near-pole azimuthal decimation
  (`Scan::Shader` has no pole-LOD path). `POLE_LOD_ENABLED` is compiled out
  under `ARDUINO`, so this capture is unaffected either way.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=SphericalHarmonics`,
`HS_PROFILE_WINDOW=16`; `just profile SphericalHarmonics` = build + flash +
capture.
