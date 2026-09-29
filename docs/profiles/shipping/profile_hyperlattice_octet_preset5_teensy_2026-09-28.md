# HyperLattice on-device profile — Octet wide (preset 5) — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot after Tiers 1 and 2 of the Octet 58 ms spec (`29305dfe2`, `f130ee632`, `1bb522211`).
Raw capture: `build/prof/octet58/preset4_ship.com4.log`. Replaces the 00:22 capture of the same name (pre-Tier-1 code). An orphaned run of the same batch took this capture from the main checkout; its provenance (source `1bb522211`, COM4, `Profile preset: 4/6`) was checked.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` base, newlib-nano, DMA LEDs, selective O3; `HyperLatticeExperimental::shade<false>` runs from cached flash as one `HS_HOT_FLASH_MEMBER` unit |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, experimental-octet-wide-flight held (`HS_PROFILE_PRESET=4`), tip `1bb522211` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 70 s capture; startup frame 1 and its window excluded; camera and palette motion live; no preset cycle or epoch crossing |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 bash tools/profile_one.sh HyperLattice profile 70 16 "-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=4"` |

Image size (instrumented single-effect image): `FLASH: code:71736, data:148748, headers:8892` / `RAM1: variables:315008, code:20168, padding:12600, free:176512` / `RAM2: variables:520064, free:4224`. The harness-attested default Phantasm image passes its size/layout gate.

Exactness cross-check: window frames 609–624 root cyc ÷ 600 MHz matches the wall sum within **0.1 ppm** (`tools/parse_profile.py ... validate`: VALID).

## Frame cadence

Frames after startup: mean render **41.686 ms**, peak **45.242 ms**, spilled **0/1096 (0.00%)** — 16 fps (every frame). Worst window (609–624) `hl_shader_draw` 42.061 ms/f.

A display window is 62.5 ms; the shader evaluates the 144×72 quadrant plus its one-pixel margin, 146×73 = 10,658 samples. `canvas_buffer_wait` is the round-up idle to the next display flip.

| | Before (spec §1) | Now |
|---|--:|--:|
| Render mean / peak | 73.299 / 76.425 ms | 41.686 / 45.242 ms |
| Shader cycles/sample (worst window) | 4,015 | 2,368 |
| Spilled | 549/549 (100%) | 0/1096 (0.00%) |

## Phase-by-phase readout

One held-preset regime.

### Held preset (window frames 609–624, worst of the capture)

```
frame                   62.490 ms  37.494 Mcyc 100.0%
  pov_preserve_half      0.138 ms   0.083 Mcyc   0.2%
  hl_shader_draw        42.061 ms  25.236 Mcyc  67.3%
  hl_timeline_step       0.013 ms   0.008 Mcyc   0.0%
  canvas_clear           0.088 ms   0.053 Mcyc   0.1%
  canvas_buffer_wait    18.427 ms  11.056 Mcyc  29.5%
```

Wall min/avg/max = 60.6/62.5/64.9 ms. The shader is the whole render; preserve, clear and timeline together take under 0.25 ms.

### Per-pixel figures

10,658 samples per frame, written directly by a premultiplied `Pixel` shader (no `filter_blend`). 2,368 shader cycles per sample in the worst window.

## Column-ISR / DMA marshaling cost

`isr_wake` 3.08%, `isr_pack` 1.57%, `isr_dma_submit` 0.21% of CPU, about 4.8% combined, leaving about 59.5 ms of foreground time per 62.5 ms window. Pack is the CPU-side LED marshaling; submit launches the asynchronous 600-byte DMA (about 400 µs on the wire at 12 MHz).

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
