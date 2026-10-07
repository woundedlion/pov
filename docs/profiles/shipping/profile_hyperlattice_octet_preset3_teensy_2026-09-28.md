# HyperLattice on-device profile — Octet 3D (preset 3) — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot after Tiers 1 and 2 of the Octet 58 ms spec (`29305dfe2`, `f130ee632`, `1bb522211`).
Raw capture: `build/prof/octet58/preset2_ship.log`. Replaces the 00:17 capture of the same name (pre-Tier-1 code, `b60fe925b`).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` base, newlib-nano, DMA LEDs, selective O3; `HyperLatticeExperimental::shade<false>` runs from cached flash as one `HS_HOT_FLASH_MEMBER` unit |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, experimental-octet-flight held (`HS_PROFILE_PRESET=2`), tip `1bb522211` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 70 s capture; startup frame 1 and its window excluded; camera and palette motion live; no preset cycle or epoch crossing |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 70 16 "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=2"` |

Image size (instrumented single-effect image): `FLASH: code:71736, data:148748, headers:8892` / `RAM1: variables:315008, code:20168, padding:12600, free:176512` / `RAM2: variables:520064, free:4224`. The harness-attested default Phantasm image passes its size/layout gate.

Exactness cross-check: window frames 913–928 root cyc ÷ 600 MHz matches the wall sum within **4.0 ppm** (`tools/parse_profile.py ... validate`: VALID).

## Frame cadence

README cells: peak 🟢 45.688, spilled 🟢 0/1096 (0.00%).

Frames after startup: mean render **39.120 ms**, peak **45.688 ms**, spilled **0/1096 (0.00%)** — 16 fps (every frame). Worst window (913–928) `hl_shader_draw` 38.772 ms/f.

A display window is 62.5 ms; the shader evaluates the 144×72 quadrant plus its one-pixel margin, 146×73 = 10,658 samples. `canvas_buffer_wait` is the round-up idle to the next display flip.

| | Before (spec §1) | Now |
|---|--:|--:|
| Render mean / peak | 69.841 / 76.566 ms | 39.120 / 45.688 ms |
| Shader cycles/sample (worst window) | 3,819 | 2,183 |
| Spilled | 549/549 (100%) | 0/1096 (0.00%) |

## Phase-by-phase readout

One held-preset regime.

### Held preset (window frames 913–928, worst of the capture)

```
frame                   62.615 ms  37.569 Mcyc 100.0%
  pov_preserve_half      0.138 ms   0.083 Mcyc   0.2%
  hl_shader_draw        38.772 ms  23.263 Mcyc  61.9%
  hl_timeline_step       0.015 ms   0.009 Mcyc   0.0%
  canvas_clear           0.088 ms   0.053 Mcyc   0.1%
  canvas_buffer_wait    21.847 ms  13.108 Mcyc  34.9%
```

Wall min/avg/max = 56.2/62.6/67.6 ms. The shader is the whole render; preserve, clear and timeline together take under 0.25 ms.

### Deep counters (`HS_PROFILE_DEEP=1`, 45 s, `build/prof/octet58/preset2_deep.log`)

Counts only; the deep image is slower and its timing is not a shipping figure.

```
frame                   62.507 ms  37.504 Mcyc 100.0%
  pov_preserve_half      0.139 ms   0.084 Mcyc   0.2%
  hl_shader_draw        44.345 ms  26.607 Mcyc  70.9%
    hl_event_step       29.975 ms  17.985 Mcyc  48.0%
      hl_layer_composite 3.274 ms   1.965 Mcyc   5.2%
      hl_event_miss      0.940 ms   0.564 Mcyc   1.5%
  canvas_buffer_wait    16.148 ms   9.689 Mcyc  25.8%
```

5.68 candidates, 4.75 misses and 0.92 layers per ray (spec §2.1 model: 5.78 / 4.86 / 0.92). `hl_event_step` counts 6.68 loop iterations per ray; one per ray is the loop exit. Per-ray setup outside the loop is about 810 cycles, and each loop iteration about 250 cycles (deep image, scope overhead included).

### Per-pixel figures

10,658 samples per frame, written directly by a premultiplied `Pixel` shader (no `filter_blend`). 2,183 shader cycles per sample in the worst window.

## Column-ISR / DMA marshaling cost

`isr_wake` 3.04% of CPU, inclusive of the nested `isr_pack` 1.54% and `isr_dma_submit` 0.21%, leaving about 60.6 ms of foreground time per 62.5 ms window before unmeasured DMA-completion and other interrupts. Pack is the CPU-side LED marshaling; submit launches the asynchronous 600-byte DMA (about 230 µs on the wire at the requested 24 MHz, LPSPI framing model).

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
