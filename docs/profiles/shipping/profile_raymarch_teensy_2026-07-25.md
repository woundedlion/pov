# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-07-25, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Raymarch`).
Raw capture: `build/prof/raymarch_ship.log`. Replaces the 2026-07-20 report;
half of a matched same-session ship/O3 pair from the full-roster 2026-07-25
sweep (tip `32478115`).

The roster index ([`../README.md`](../README.md)) ranks this row from a
2026-07-26 11:38 log that postdates this report: the index peak is the current
figure and the numbers below are the earlier capture. `just profile Raymarch`
regenerates this report against a fresh log.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` baseline + `HS_O3` regions. The volume march path (`Scan::Volume`, inside `rm_shader_draw`) is in an HS_O3 region, including the `always_inline`'d `quintic_kernel`. |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Raymarch 288×144, single-entry playlist, tip `32478115` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Raymarch profile 70 32` |

Image size: `FLASH: code:94900, data:1068796, headers:8776` / `RAM1:
variables:314368, code:38520, padding:27016, free:144384` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 161–192 root counter cyc ÷ 600 MHz matches
the measured wall sum within **3.1 ppm**.

## Frame cadence

**Pass aggregate**: peak frame render **53.02 ms** (frames 1–32), spilled
**0/1088 frames (0%)** — 🟢. Render averages ~47 ms with peaks near 53 ms;
`rm_buffer_wait` ~24% — holds 16 fps with ~9 ms of headroom.

## Phase-by-phase readout

Phase schedule: Raymarch continuously marches a volume SDF; single regime.

### Peak window (frames 1–32)

```
frame              61.68 ms  37.01 Mcyc  100%
  rm_timeline_step 46.84 ms  28.11 Mcyc   75%
    rm_shader_draw 46.80 ms  28.08 Mcyc   99%
      filter_blend  0.87 ms   0.53 Mcyc    1%  x5548  95 cyc/blend
  rm_buffer_wait   14.84 ms   8.90 Mcyc   24%
```

Wall min/avg/max = 49.05/61.68/71.85 ms; percentages are of the parent scope.
`rm_shader_draw` is the whole render — the per-pixel volume march — with a
negligible blend tail (5,548 blended px = 0.54× the quadrant).

### Per-pixel figures

Quadrant = 10,368 px; `filter_blend` = 5,548 blended px/frame ⇒ **0.54×
coverage**. `rm_shader_draw` 28.08 Mcyc ÷ 10,368 = **2,708 cyc per pixel** (the
march is evaluated per canvas pixel, not per blended pixel).

## Column-ISR / DMA marshaling cost

```
isr_wake        1397/frame  min/avg/max 0.51/1.53/20.91 us  cpu 2.81%
isr_pack         139/frame  min/avg/max 6.11/6.57/9.42 us  cpu 1.20%
isr_dma_submit   139/frame  min/avg/max 0.60/0.93/1.09 us  cpu 0.16%
```

(window frames 1–32.) Total ISR CPU share **4.17%** — the lowest on the roster,
because Raymarch's longer frames mean fewer packs/submits per unit render.

## Summary ranking

1. `rm_shader_draw` — 99% of `rm_timeline_step`, 46.8 ms: the per-pixel volume
   march. Fully inside the `Scan::Volume` HS_O3 region.

Raymarch holds a comfortable 16 fps (~9 ms margin). The 2026-07-20 landing
(66.96 → 53.34 ms) removed divides/sqrts from the near-field SDF and added an
`always_inline` on `quintic_kernel`; this capture reproduces that at 53.02 ms.
The global-`-O3` ceiling buys only 53.02 → 52.76 ms (1.005×) precisely because
the march path is **already** `-O3` in the shipping image.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `rm_shader_draw`; coverage is sparse.
- Selective -O3: the volume march is the archetype covered region — which is why
  the shipping and ceiling images nearly coincide.
- Working tree tip `32478115`, only untracked docs.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Raymarch`,
`HS_PROFILE_WINDOW=32`; `just profile Raymarch` = build + flash + capture.
