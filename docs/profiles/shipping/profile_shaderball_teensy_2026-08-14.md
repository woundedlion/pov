# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-14, **selective -O3**)

Final 59 ms acceptance capture of the current 12-preset shipping bank. Raw
captures:
`build/prof/shaderball_shipping_full_c06c5fbc_com3.log` and
`build/prof/shaderball_shipping_full_c06c5fbc_com4.log`. Their adjacent
`.provenance` files attest the source, compiler, ELF hashes, and archived
artifacts.

## Setup

| | |
|---|---|
| Hardware | Two Teensy 4.0 boards at 600 MHz, segmented POV driver, COM3 and COM4 |
| Image | `profile`: `-Os` base plus shipping selective-O3 regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, source `c06c5fbc` |
| Compiler | Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1 20251203 |
| Method | Two independent 650 s captures, 16-frame windows, normal 480-revolution choreography |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 650 16` with `HS_TEENSY_PORT` pinned per board |

Both captures contain 648 windows and 10,368 frames. Each visits presets 0–11,
wraps 11→0, has monotonic frame numbers, no epoch reset, and zero spilled
frames. Root cycles divided by 600 MHz agree with wall time within 3.5 ppm on
COM3 and 3.7 ppm on COM4.

The profile image is:

```text
FLASH  code 113,656 B  data 168,900 B  headers 8,260 B
RAM1   variables 323,200 B  code 28,008 B  padding 4,760 B
RAM1   local free 168,320 B
RAM2   variables 520,064 B  allocator free 4,224 B
```

COM3 profile ELF SHA-256 is
`85587dafd3300d7e1c59303c2d6cbeeefe311c3dad644fc8d5cfb6874c423ec3`;
COM4 is
`dae7569d441b1f49ea95b19a2117995d23676a1163e85e778ff3aca2108b6c49`.

## Result

Every current preset and transition is below **59.00 ms** on both boards.
COM3 and COM4 both peak at **55.69 ms** in preset 0, leaving 3.31 ms of hard
gate headroom. Each capture reports **0/10,368 spilled frames**.

| Preset | Authored topology | COM3 peak, ms | COM4 peak, ms | Worst, ms | Clean shader, ms |
|---:|---|---:|---:|---:|---:|
| 0 | Glitch wave-shear grid | 55.69 | 55.69 | 55.69 | 48.33 |
| 11 | Gnomonic dodecahedral vector-noise + mirror grid | 54.76 | 54.76 | 54.76 | 46.95 |
| 6 | Gnomonic dodecahedral wave/mirror grid | 52.23 | 52.23 | 52.23 | 45.11 |
| 5 | Peirce dodecahedral grid | 51.84 | 51.74 | 51.84 | 44.35 |
| 10 | Stereographic prism polar-wave lattice | 48.28 | 48.23 | 48.28 | 40.57 |
| 9 | Sinusoidal curl lattice, coarse field | 48.04 | 47.98 | 48.04 | 40.58 |
| 8 | Sinusoidal curl lattice, fine field | 47.58 | 47.60 | 47.60 | 40.55 |
| 7 | Gnomonic affine lattice contour | 44.94 | 44.94 | 44.94 | 38.02 |
| 4 | Bonne kaleidoscope lattice + mirror | 41.30 | 41.34 | 41.34 | 38.44 |
| 3 | Gnomonic glitch grid + mirror | 25.86 | 25.87 | 25.87 | 22.69 |
| 1 | Kaleidoscope twin-wave + inner mirror | 25.76 | 25.72 | 25.76 | 22.70 |
| 2 | Gnomonic kaleidoscope grid + mirror | 25.36 | 25.39 | 25.39 | 22.71 |

The peak includes normal through-clear transitions. The clean-shader column is
the worse board's maximum clean-hold `sb_shader_draw` window.

## Worst-window readout

COM3 frames 5361–5376 contain the pass peak and a full preset-0 hold:

```text
frame                 62.39 ms  37.43 Mcyc  100%
  sb_shader_draw      48.33 ms  29.00 Mcyc   77%
  canvas_buffer_wait   7.37 ms   4.42 Mcyc   11%
  canvas_clear        92.6 us   55.6 kcyc     0%
  sb_timeline_step    69.5 us   41.7 kcyc     0%
```

Wall min/avg/max is 60.85/62.39/63.85 ms; render averages 55.02 ms and peaks
at 55.69 ms. The clean shader consumes about 2,797 cycles per visible sample
across the 10,368-sample quadrant.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.64/1.96/13.78 us  cpu 3.62%
isr_pack         144/frame  min/avg/max 6.29/7.10/9.38 us   cpu 1.63%
isr_dma_submit   144/frame  min/avg/max 0.61/0.94/6.11 us   cpu 0.21%
```

ISR work is included in the cycle counters. Pack plus submit costs about
1.16 ms of CPU per rendered frame in this window; LED transfer continues
asynchronously after submission.

## Historical baseline and roster remap

The 2026-08-13 17-preset baseline peaked at 147.35 ms with 664/1,472 spills.
The optimized pre-reduction capture at `c5cb0bb4` peaked at 58.23 ms with no
spills. Those measurements keep their original indices and are not relabelled
as current-preset results.

| Former preset | Current preset |
|---:|---:|
| 0 | 0 |
| 1 | 1 |
| 2 | 2 |
| 3 | 3 |
| 4 | 4 |
| 5 | Retired |
| 6 | Retired |
| 7 | Retired |
| 8 | 5 |
| 9 | Retired |
| 10 | Retired |
| 11 | 6 |
| 12 | 7 |
| 13 | 8 |
| 14 | 9 |
| 15 | 10 |
| 16 | 11 |

The current peak is 3.31 ms below the gate and 2.54 ms below the historical
pre-reduction peak. The five retired programs also remove the former 0.77 ms
flash-layout sentinel from the shipping bank.

## Matched global-O3 comparison

The [source-matched global-O3 report](../O3/profile_shaderball_teensy_2026-08-14.md)
uses the same source and two-board protocol. Global O3 peaks at 84.64 ms,
spills 1,440/8,928 frames on each board, and has five red presets. It is
rejected for ShaderBall; selective O3 is both the shipping configuration and
the faster live result.

The deduplicated ShaderBall-named symbol audit totals 41,076 B in shipping
(2,806 B ITCM and 38,270 B flash) versus 43,224 B under global O3 (3,782 B
ITCM and 39,442 B flash). Global O3 adds 2,148 B of named code while regressing
the device peak by 28.95 ms. Whole profile-image deltas are +13,000 B flash
code and +12,352 B ITCM.

## Full-roster size and validation

The current Phantasm image passes every memory gate:

```text
FLASH code          422,560 B  (-5,416 B vs pre-reduction)
RAM1 variables      314,880 B
ITCM code           195,144 / 196,608 B  (1,464 B remaining)
RAM1 padding          1,464 B
RAM1 local free      12,800 B
RAM2 variables      520,064 B
RAM2 allocator free   4,224 B
```

The roster reduction reclaimed 1,296 B of ITCM and 5,416 B of flash. RAM1 and
RAM2 floors are unchanged. The COM3 preflight Phantasm ELF SHA-256 is
`3be67928a826425f1ecc9bd9bcea9b17f756b4080b160b02b797661ec62daff1`;
the complete provenance and artifact paths are in
`build/prof/shaderball_shipping_full_c06c5fbc_com3.provenance`.

## Caveats

- Counter scopes include live flywheel/DMA ISR time.
- Preset buckets own their on-screen frames, including transition work.
- The runs use normal choreography, not an accelerated diagnostic cycle.
- Compiler, linker, manifest, or flash-layout changes require another
  two-board full cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShaderBall`, and
`HS_PROFILE_WINDOW=16`; `tools/profile_one.sh` performs locked build, flash,
capture, marker checks, and artifact attestation.
