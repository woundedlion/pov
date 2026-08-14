# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-14, **-O3**)

Source-matched global-O3 reference for the current 12-preset bank. The
[shipping selective-O3 report](../shipping/profile_shaderball_teensy_2026-08-14.md)
is the acceptance result. Raw captures:
`build/prof/shaderball_o3_full_c06c5fbc_com3.log` and
`build/prof/shaderball_o3_full_c06c5fbc_com4.log`. Their adjacent
`.provenance` files attest the source, compiler, ELF hashes, and archived
artifacts.

## Setup

| | |
|---|---|
| Hardware | Two Teensy 4.0 boards at 600 MHz, segmented POV driver, COM3 and COM4 |
| Image | `profile_o3`: global `-O3 -ffast-math` reference image |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, source `c06c5fbc` |
| Compiler | Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 20251203 |
| Method | Two independent 650 s captures, 16-frame windows, normal 480-revolution choreography |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile_o3 650 16` with `HS_TEENSY_PORT` pinned per board |

Both captures contain 558 windows and 8,928 frames. Each visits presets 0–11,
wraps 11→0, has monotonic frame numbers, and has no epoch reset. Root cycles
divided by 600 MHz agree with wall time within 0.4 ppm on COM3 and 2.1 ppm on
COM4.

The profile image is:

```text
FLASH  code 126,656 B  data 169,240 B  headers 8,232 B
RAM1   variables 323,200 B  code 40,360 B  padding 25,176 B
RAM1   local free 135,552 B
RAM2   variables 520,064 B  allocator free 4,224 B
```

COM3 profile ELF SHA-256 is
`8fe303b25ac0f69173d7802d219f559ee2fd347fd0a0a9b8281b3df37e390805`;
COM4 is
`832ec8d42360619aba7c6a3f1ecc99d60d1d1ce3f0ff744b9cd755d502d882e2`.

## Result

Global O3 fails the 62.5 ms display-window gate on five presets. COM3 peaks at
**84.64 ms** and COM4 at **84.39 ms**, both in preset 11. Each capture spills
**1,440/8,928 frames (16.1%)**. The five red presets account for 1,440/2,405
of their own frames (59.9%).

| Preset | Authored topology | COM3 peak, ms | COM4 peak, ms | Worst, ms | Clean shader, ms | Cadence |
|---:|---|---:|---:|---:|---:|---|
| 11 | Gnomonic dodecahedral vector-noise + mirror grid | 84.64 | 84.39 | 84.64 | 76.51 | Red |
| 0 | Glitch wave-shear grid | 83.95 | 84.20 | 84.20 | 76.70 | Red |
| 9 | Sinusoidal curl lattice, coarse field | 69.28 | 68.94 | 69.28 | 61.22 | Red |
| 10 | Stereographic prism polar-wave lattice | 68.39 | 68.28 | 68.39 | 60.93 | Red |
| 8 | Sinusoidal curl lattice, fine field | 65.57 | 65.79 | 65.79 | 57.75 | Red |
| 2 | Gnomonic kaleidoscope grid + mirror | 59.98 | 60.54 | 60.54 | 55.18 | Green |
| 3 | Gnomonic glitch grid + mirror | 59.16 | 59.74 | 59.74 | 54.72 | Green |
| 6 | Gnomonic dodecahedral wave/mirror grid | 58.89 | 58.94 | 58.94 | 50.99 | Green |
| 7 | Gnomonic affine lattice contour | 58.67 | 58.82 | 58.82 | 50.42 | Green |
| 5 | Peirce dodecahedral grid | 49.56 | 49.62 | 49.62 | 41.89 | Green |
| 4 | Bonne kaleidoscope lattice + mirror | 41.03 | 41.00 | 41.03 | 37.27 | Green |
| 1 | Kaleidoscope twin-wave + inner mirror | 31.44 | 31.45 | 31.45 | 27.91 | Green |

The peak includes normal transition work. The clean-shader column is the worse
board's maximum clean-hold `sb_shader_draw` window.

## Worst-hold readout

COM3 frames 5057–5072 are a full preset-11 hold adjacent to the pass peak:

```text
frame                124.78 ms  74.87 Mcyc  100%
  sb_shader_draw      76.50 ms  45.90 Mcyc   61%
  canvas_buffer_wait  41.20 ms  24.72 Mcyc   33%
  canvas_clear        89.7 us   53.8 kcyc     0%
  sb_timeline_step    63.9 us   38.4 kcyc     0%
```

Wall min/avg/max is 123.86/124.78/125.40 ms; render averages 83.58 ms and
peaks at 84.20 ms in this hold. The pass peak of 84.64 ms occurs in the
preceding transition window. The clean shader consumes about 4,428 cycles per
visible sample.

## Column-ISR / DMA marshaling cost

```text
isr_wake        2302/frame  min/avg/max 0.55/1.82/21.74 us  cpu 3.35%
isr_pack         288/frame  min/avg/max 6.32/6.87/9.40 us   cpu 1.58%
isr_dma_submit   288/frame  min/avg/max 0.61/0.91/9.90 us   cpu 0.21%
```

The two-window cadence doubles ISR calls per rendered frame. Pack plus submit
costs about 2.24 ms of CPU per frame; LED transfer remains asynchronous.

## Global O3 versus selective O3

The matched shipping image peaks at 55.69 ms with zero spills on both boards.
Global O3 regresses the worst frame by **28.95 ms (52.0%)** and turns five
presets red. It is rejected as a ShaderBall optimization and retained only as
a compiler reference.

The deduplicated ShaderBall-named symbol audit is:

| Image | ITCM | Flash | Total |
|---|---:|---:|---:|
| Shipping selective O3 | 2,806 B | 38,270 B | 41,076 B |
| Global O3 | 3,782 B | 39,442 B | 43,224 B |
| O3 minus shipping | +976 B | +1,172 B | +2,148 B |

The whole profile image adds 13,000 B of flash code and 12,352 B of ITCM
under global O3. The full-roster Phantasm gate remains a shipping selective-O3
measurement; the profile image is not a memory-acceptance proxy.

## Caveats

- Counter scopes include live flywheel/DMA ISR time.
- Preset buckets own their on-screen frames, including transition work.
- The runs use normal choreography, not an accelerated diagnostic cycle.
- This global-O3 image is a reference build, not a shipping candidate.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShaderBall`, and
`HS_PROFILE_WINDOW=16`; use
`bash tools/profile_one.sh ShaderBall profile_o3 650 16` for this reference
configuration.
