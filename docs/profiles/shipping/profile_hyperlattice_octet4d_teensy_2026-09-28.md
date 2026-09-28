# HyperLattice Octet 4D slice (preset 4) on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot after Tiers 1 and 2 of the Octet 58 ms spec (`29305dfe2`, `f130ee632`, `1bb522211`).
Raw capture: `build/prof/octet58/preset3_ship.com3.log`. Supersedes the 4D rows of [the 2026-09-27 octet report](profile_hyperlattice_octet_teensy_2026-09-27.md), which stays as the pre-optimization record the optimization ledger cites. Captured on COM3 (identical board) from a worktree at the same commit.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` base, newlib-nano, DMA LEDs, selective O3; `HyperLatticeExperimental::shade<true>` runs from cached flash as one `HS_HOT_FLASH_MEMBER` unit |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, experimental-octet-4d-slice held (`HS_PROFILE_PRESET=3`), tip `1bb522211` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 45 s capture; startup frame 1 and its window excluded; camera and palette motion live; no preset cycle or epoch crossing |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice profile 45 16 "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=3"` |

Image size (instrumented single-effect image): `FLASH: code:71736, data:148748, headers:8892` / `RAM1: variables:315008, code:20168, padding:12600, free:176512` / `RAM2: variables:520064, free:4224`. The harness-attested default Phantasm image passes its size/layout gate.

Exactness cross-check: window frames 145–160 root cyc ÷ 600 MHz matches the wall sum within **0.3 ppm** (`tools/parse_profile.py ... validate`: VALID).

## Frame cadence

Frames after startup: mean render **137.407 ms**, peak **143.318 ms**, spilled **233/233 (100%)** — 5.33 fps (3 display windows per frame). Worst window (145–160) `hl_shader_draw` 137.256 ms/f.

A display window is 62.5 ms; the shader evaluates the 144×72 quadrant plus its one-pixel margin, 146×73 = 10,658 samples. `canvas_buffer_wait` is the round-up idle to the next display flip.

| | Before (spec §1) | Now |
|---|--:|--:|
| Render mean / peak | 177.308 / 188.387 ms | 137.407 / 143.318 ms |
| Shader cycles/sample (worst window) | 9,865 | 7,727 |
| Spilled | 231/231 (100%) | 233/233 (100%) |

## Phase-by-phase readout

One held-preset regime.

### Held preset (window frames 145–160, worst of the capture)

```
frame                  187.504 ms 112.502 Mcyc 100.0%
  pov_preserve_half      0.133 ms   0.080 Mcyc   0.1%
  hl_shader_draw       137.256 ms  82.353 Mcyc  73.2%
  hl_timeline_step       0.006 ms   0.004 Mcyc   0.0%
  canvas_clear           0.086 ms   0.051 Mcyc   0.0%
  canvas_buffer_wait    48.263 ms  28.958 Mcyc  25.7%
```

Wall min/avg/max = 186.8/187.5/188.1 ms. The shader is the whole render; preserve, clear and timeline together take under 0.25 ms.

### Deep counters (`HS_PROFILE_DEEP=1`, 45 s, `build/prof/octet58/preset3_deep.log`)

Counts only; the deep image is slower and its timing is not a shipping figure.

```
frame                  187.455 ms 112.473 Mcyc 100.0%
  pov_preserve_half      0.131 ms   0.079 Mcyc   0.1%
  hl_shader_draw       143.874 ms  86.324 Mcyc  76.8%
    hl_event_step      105.306 ms  63.183 Mcyc  56.2%
      hl_layer_composite 2.228 ms   1.337 Mcyc   1.2%
      hl_event_miss      2.371 ms   1.423 Mcyc   1.3%
  canvas_buffer_wait    41.564 ms  24.939 Mcyc  22.2%
```

10.08 candidates, 9.25 misses and 0.84 layers per ray (spec §4.1 model with ownership: 10.06 candidates; 14.11 before). Ownership removed the candidates the model predicted, but per-candidate cost did not fall to the modelled twenty or so instructions per class: each loop iteration costs about 535 cycles, and per-ray setup outside the loop (ownership for 12 classes and the 4D projection) about 2,170 cycles, 27% of the shader.

### Per-pixel figures

10,658 samples per frame, written directly by a premultiplied `Pixel` shader (no `filter_blend`). 7,727 shader cycles per sample in the worst window.

## Column-ISR / DMA marshaling cost

`isr_wake` 3.04%, `isr_pack` 1.55%, `isr_dma_submit` 0.21% of CPU, about 4.8% combined, leaving about 59.5 ms of foreground time per 62.5 ms window. Pack is the CPU-side LED marshaling; submit launches the asynchronous 600-byte DMA (about 400 µs on the wire at 12 MHz).

## Summary ranking

1. `hl_shader_draw`: plane-crossing traversal, strut coverage and compositing; the entire render.
2. `canvas_buffer_wait`: display alignment idle.

## Caveats

- CYCCNT includes ISR time. Startup frame 1 and its window are excluded.
- Fixed preset pinning keeps camera and palette motion live; results cover the recorded trajectory.
- No global-O3 twin was captured for this code.
- The experimental presets are opt-in (`HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1`).

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16` and `HS_PROFILE_PRESET=<index>`; command in Setup.
