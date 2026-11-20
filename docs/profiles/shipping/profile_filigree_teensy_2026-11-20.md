# Filigree on-device profile — Teensy 4.0, segmented mode (2026-11-20, **selective -O3**)

Raw capture `C:/work/temp/holosphere-20261119-filigree-profile/morph-shared-cell/filigree_ship.log`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0 on COM3, 600 MHz |
| Image | `profile`, shipping `-Os` with selective `HS_O3` |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master with live flywheel and DMA ISRs |
| Effect | Filigree, 288×144, source `c98ca665517a71b480f8fc6f4af1ec907681a724` |
| Method | 176 s capture, 16 frames per readout |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh Filigree profile 176 16` |

Captured 2026-11-20 08:42, America/Los_Angeles. Compiler `GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`.
Raw SHA-256 `fd46bb9db1eb3dbc4cc6e75ec1f8de7b571995c99f57eba44b32fb4de3de1073`.
Profile ELF SHA-256 `7c031fccb372ed3702af8ca024d0aaa595dedafdde816321a4c8119bb5db0dc9`.
Phantasm ELF SHA-256 `c3621b48a9c5ce5350c86f94001dfab4eb478e8798b7064a96a93823fa5b66fe`.

Profile image size `FLASH: code:116200, data:456380, headers:9052   free for files:1449984 / RAM1: variables:314144, code:34424, padding:31112   free for local variables:144608 / RAM2: variables:520064  free for malloc/new:4224`.

Shipping Phantasm image size `FLASH: code:1115344, data:724920, headers:9072   free for files:182280 / RAM1: variables:313760, code:195832, padding:776   free for local variables:13920 / RAM2: variables:520064  free for malloc/new:4224`.
Compared with the previous combined Filigree build, FLASH code changed by
+4,448 B and FLASH data by +107,096 B. RAM1 code changed by
+432 B and variables by +0 B.

## Frame cadence

Setup frame 1 took **76.754 ms** to render. Both steady-state holds and
both transition directions passed the 62.5 ms display budget. Across
**2,791 live frames**, the peak was **60.862 ms**
on frame 2159 (Ember Lattice → Jade Weave). Peak foreground wall time was
60.862 ms. Every live render and foreground wall measurement
stayed within budget.

The review extended beyond Filigree's 105 s playlist slot to include the
complete reverse transition. Within the scheduled slot, **0/1,667**
live frames exceeded budget. Peak render was **57.375 ms**
(Jade Weave → Ember Lattice).

README cells: peak 🟢 57.375 (2), spilled 🟢 0/1667 (0.00%).

README setup frame 1: 76.754 ms.

## Phase-by-phase readout

Each authored dwell is 600 frames. Each transition takes 480 sinusoidally
eased steps. Each transition range includes the incoming endpoint frame.

| Phase | Frames | Count | Mean render ms | Peak render ms | Peak frame | Over budget |
|---|---:|---:|---:|---:|---:|---:|
| Jade Weave hold | 2–599 | 598 | 41.152 | 42.001 | 94 | 0 |
| Jade Weave → Ember Lattice | 600–1080 | 481 | 54.446 | 57.375 | 1080 | 0 |
| Ember Lattice hold | 1081–1678 | 598 | 34.246 | 35.305 | 1401 | 0 |
| Ember Lattice → Jade Weave | 1679–2159 | 481 | 53.755 | 60.862 | 2159 | 0 |
| Jade Weave hold | 2160–2757 | 598 | 41.282 | 42.683 | 2561 | 0 |
| Jade Weave → Ember Lattice | 2758–2792 | 35 | 53.844 | 55.852 | 2785 | 0 |

The previous combined effect peaked at 71.527 ms with a discrete solid switch
on frame 1080. The first source-morph build peaked at 70.021 ms with two
misses at the final transition frames.

Filigree blends both edge fields in the source stage. Packed narrow lookups
and bounded short-angle kernels preserve oriented arc coordinates and coverage.
Both solids share one cubemap cell address, and certified single-edge interiors
bypass the candidate search. The scan factors the camera rotation per column and skips color work for
coverage that rounds to black in RGB16. The phase stack samples and blends
both palette endpoints directly into a single color carrier, preserving
RGB-first interpolation before perceptual modulation. Jade Weave's return
frame defers the new hue table while neither live layer rotates hue.

## Summary ranking

Frames 2145–2160. A clipped frame shades 10,368
pixels. Render includes interrupt time and excludes the buffer wait.

| Counter | ms/frame | Share of foreground | Calls/frame |
|---|---:|---:|---:|
| `fx_shader_draw` | 52.979 | 98.4% | 1.00 |
| `fx_prepare_column` | 5.617 | 10.4% | 144.00 |
| `fx_palette_base_bake` | 0.167 | 0.3% | 0.94 |
| `pov_preserve_half` | 0.132 | 0.2% | 1.00 |
| `fx_prepare_frame` | 0.115 | 0.2% | 1.00 |
| `fx_timeline_step` | 0.108 | 0.2% | 1.00 |
| `canvas_clear` | 0.095 | 0.2% | 1.00 |
| `fx_advance` | 0.077 | 0.1% | 1.00 |
| `fx_palette_rgb_bake` | 0.076 | 0.1% | 0.06 |

### Column ISR and DMA

| Counter | Calls/frame | Min/mean/max μs | CPU |
|---|---:|---:|---:|
| `isr_wake` | 1138.1 | 0.505/1.668/26.991 | 3.07% |
| `isr_pack` | 142.2 | 5.406/6.656/9.551 | 1.53% |
| `isr_dma_submit` | 142.2 | 0.628/0.915/6.340 | 0.21% |

The capture used the shipping preset order, dwell and easing. The stock
harness measured the full-width setup frame separately. The capture used
the segmented hardware driver with active ISR and DMA and covered Filigree's
standalone render windows.

## Harness

The stock profile comes from `targets/Profile/Profile.ino` with
`HS_PROFILE_TARGET=Filigree` and `HS_PROFILE_WINDOW=16`. The command above
acquires the board through `tools/profile_one.sh`.
