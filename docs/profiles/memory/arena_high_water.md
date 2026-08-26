# Arena high-water per effect (host probe)

Measured 2026-07-15 with `tests/arena_measure.cpp` (host Clang `-Os`,
288×144, `init()` + 8 frames per effect, 8 MiB host arena). Values are
`Arena::get_high_water_mark()` for the three partitions of the one global
block. Two snapshots below: the survey that motivated the trail
compression, and the post-landing state. Regenerate with:

```
cmake --build build/tests --target arena_measure && build/tests/tests/arena_measure.exe
```

Host figures are a conservative **upper bound** on the device: the host is
64-bit, so pointer-bearing pooled headers inflate vs the 32-bit device (see
the arena_measure.cpp header). Two further caveats: the probe samples one
random spawn per effect (shape/preset draws vary run to run), and 8 frames
does not capture long-run growth — this is the same window the CI arena
budget gate trusts for device sizing.

Each snapshot lists the roster of its own commit. Flyby and Liquid2D below
predate their merge into ShaderWorkbench.

## Survey snapshot (pre trail-compression, master 2b875b6e)

| Effect | persist | scratch A | scratch B | total |
|---|--:|--:|--:|--:|
| MindSplatter | 327,120 | 3,684 | 0 | 330,804 |
| GSReactionDiffusion | 175,104 | 122,880 | 0 | 297,984 |
| BZReactionDiffusion | 164,352 | 115,200 | 0 | 279,552 |
| DisplacementField | 242,404 | 0 | 0 | 242,404 |
| Dynamo | 130,880 | 2,496 | 0 | 133,376 |
| Fishbowl | 30,284 | 96,528 | 0 | 126,812 |
| HopfFibration | 108,072 | 2,688 | 0 | 110,760 |
| IslamicStars | 29,696 | 18,376 | 13,438 | 61,510 |
| RingShower | 49,152 | 0 | 0 | 49,152 |
| Comets | 33,256 | 0 | 1,608 | 34,864 |
| GnomonicStars | 27,112 | 0 | 1,608 | 28,720 |
| Voronoi | 14,400 | 13,750 | 0 | 28,150 |
| HankinSolids | 16,216 | 8,208 | 288 | 24,712 |
| DreamBalls | 14,568 | 6,320 | 2,798 | 23,686 |
| RingSpin | 18,176 | 0 | 1,608 | 19,784 |
| Thrusters | 0 | 16,224 | 0 | 16,224 |
| MeshFeedback | 10,528 | 3,072 | 1,016 | 14,616 |
| Raymarch | 3,072 | 3,340 | 4,138 | 10,550 |
| PetalFlow | 0 | 9,216 | 0 | 9,216 |
| MobiusGrid | 3,112 | 5,856 | 0 | 8,968 |
| ShapeShifter | 3,072 | 2,670 | 1,608 | 7,350 |
| SphericalHarmonics | 3,072 | 0 | 1,608 | 4,680 |
| Flyby | 3,072 | 0 | 0 | 3,072 |
| Liquid2D | 3,072 | 0 | 0 | 3,072 |

Per-arena maxima: persist 327,120 · scratch A 122,880 · scratch B 13,438.
Worst single-effect total: **330,804 B (MindSplatter)** against the
`DEVICE_GLOBAL_ARENA_SIZE` budget of 337,920 B — 97.9 % of the budget, a
7,116 B margin.

## Post-landing state (master b9c5fa6a, snorm16 trail compression)

`animation: quantize Particle trail history to snorm16` shrank
MindSplatter's pool from 316 to 180 B/particle. Only its row changed:

| Effect | persist | scratch A | scratch B | total |
|---|--:|--:|--:|--:|
| MindSplatter | 187,856 | 3,684 | 0 | 191,540 |

Worst single-effect total is now **297,984 B (GSReactionDiffusion)**.

## FlexRAM bank: freed (master aeba37b5)

`memory: shrink the device arena to 298 KiB` set
`DEVICE_GLOBAL_ARENA_SIZE` to 305,152 B — GSReactionDiffusion fits with
7,168 B margin, every other effect with ≥ 25 KB. Measured on the phantasm
ELF: RAM1 variables 312,704 B, free-for-locals **47,744 B** (was 14,976) —
one whole 32 KiB bank of headroom for ITCM -O3 promotions. The size-gate
arena bracket moved down one bank to [294,912, 327,680].

On-device MindSplatter profile (2026-07-16, full 5-preset cycle, both
configs): the snorm16 trail decode is timing-neutral in both images —
shipping preset-0 draw scope 96.49 → 97.30 ms (+0.84 %) vs the 2026-07-14
baseline, and the day-apart global-O3 A/B (pre- vs post-compression)
reproduces within ~2 % per preset, inside the ±10 % per-preset coverage
noise. Cadence class unchanged (8 fps steady; worst preset 108.57 ms ship /
88.92 ms O3). That capture has since been superseded; the current reports are
[shipping](../shipping/profile_mindsplatter_teensy_2026-08-26.md) and
[O3](../O3/profile_mindsplatter_teensy_2026-08-26.md).

A second bank (budget 272,384 B) would additionally require shrinking both
reaction-diffusion effects, whose footprints are mutable simulation grids —
not realistic via data-format changes alone.
