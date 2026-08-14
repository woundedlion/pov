# ShaderBall on-device profile — Teensy 4.0, segmented mode (2026-08-14, selective O3)

Final 59 ms quest capture of the curated 17-preset shipping bank. Raw captures:
`build/prof/shaderball_shipping_full_c5cb0bb4_com3.log` and
`build/prof/shaderball_shipping_full_c5cb0bb4_com4.log`.

## Setup

| | |
|---|---|
| Hardware | Two Teensy 4.0 boards at 600 MHz, segmented POV driver, COM3 and COM4 |
| Image | `profile`: `-Os` base plus shipping selective-O3 regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShaderBall 288×144, single-entry playlist, source `c5cb0bb4` |
| Method | Two independent 650 s captures, 16-frame windows, normal 480-revolution choreography |
| Reproduce | `bash tools/profile_one.sh ShaderBall profile 650 16` with `HS_TEENSY_PORT` pinned per board |

Both 648-window captures have monotonic frame numbers, visit all 17 presets,
wrap 16→0, and contain no epoch reset. Root cycles divided by 600 MHz agree
with wall time within 2.8 ppm.

## Result

Every shipping preset and transition is below **59.00 ms** on both boards.
COM3 peaks at **58.23 ms** and COM4 at **58.06 ms**, both in preset 9, with
**0/10,368 spilled frames** in each capture. The 2026-08-13 baseline peaked at
147.35 ms and spilled 664/1,472 frames; the final peak is 60.5% lower.

| Preset | Authored topology | COM3 peak, ms | COM4 peak, ms | Worst, ms |
|---:|---|---:|---:|---:|
| 9 | Dodecahedral noise grid | 58.23 | 58.06 | 58.23 |
| 8 | Peirce dodecahedral grid | 57.75 | 57.78 | 57.78 |
| 16 | Gnomonic dodecahedral vector-noise + mirror grid | 56.19 | 56.15 | 56.19 |
| 0 | Glitch wave-shear grid | 55.43 | 55.42 | 55.43 |
| 11 | Gnomonic dodecahedral wave/mirror grid | 53.77 | 53.75 | 53.77 |
| 10 | Dodecahedral noise lattice + mirror | 53.59 | 53.67 | 53.67 |
| 7 | Dodecahedral grid + mirror | 51.15 | 51.20 | 51.20 |
| 14 | Sinusoidal curl lattice, coarse field | 48.94 | 49.04 | 49.04 |
| 15 | Stereographic prism polar-wave lattice | 48.89 | 48.88 | 48.89 |
| 13 | Sinusoidal curl lattice, fine field | 48.39 | 48.41 | 48.41 |
| 12 | Gnomonic affine lattice contour | 45.56 | 45.56 | 45.56 |
| 6 | Kaleidoscope noise grid + edge fade | 45.33 | 45.38 | 45.38 |
| 5 | Peirce kaleidoscope lattice | 42.28 | 42.31 | 42.31 |
| 4 | Bonne kaleidoscope lattice + mirror | 41.59 | 41.64 | 41.64 |
| 1 | Kaleidoscope twin-wave + inner mirror | 26.87 | 26.84 | 26.87 |
| 3 | Gnomonic glitch grid + mirror | 25.90 | 26.07 | 26.07 |
| 2 | Gnomonic kaleidoscope grid + mirror | 25.53 | 25.53 | 25.53 |

Preset 9 has only 0.77 ms of hard-gate headroom and preset 8 has 1.22 ms.
Their independent two-board repeats pass, but they remain the first regression
sentinels after compiler or flash-layout changes.

## What closed the quest

The retained implementation combines several exact structural reductions:

- analytic Simplex curl derivatives, a specialized curl surface path, and a
  fused Euler step;
- coupled three-channel vector noise and prepared vector-warp trigonometry;
- prepared wave-shear direction and a persistent spherical hue-noise field;
- removal of redundant palette lookup and hidden palette rebakes;
- visible-endpoint-only palette work and stopping discarded transition runtime;
- selective ITCM residence for presets 0, 8, 10, and 11, funded by moving
  cold control/setup code to flash; and
- 4 KiB-aligned flash wrappers for the phase-sensitive preset-9 and preset-4
  kernels.

The complete per-experiment ledger, including rejected and superseded A/Bs,
uses the required `SHA / Experiment / Preset/run / Peak before -> after /
Ship/O3 symbol delta / Phantasm ITCM/RAM delta / Tests / Decision` schema in
the [quest plan](../../shaderball_59ms_quest_plan.md#quest-outcome-and-experiment-ledger).

The last item was established by generated-code A/B rather than source-level
guessing. The final Phantasm ELF places the 0x37c-byte preset-9 shader at
`0x60035000` and the 0x300-byte preset-4 shader at `0x60036000`. The resident
shader bodies occupy ITCM at `0x5a4` (preset 11, 0x2d4 bytes), `0x878`
(preset 8, 0x2b4), `0xb2c` (preset 0, 0x298), and `0xdc4` (preset 10, 0x42c).

## Full-roster size and validation

The attested Phantasm image passes every memory gate:

```text
FLASH code       427,976 B
RAM1 variables   314,880 B
ITCM code        196,440 / 196,608 B (168 B remaining)
RAM1 local free   12,800 B
RAM2 allocator free 4,224 B
```

Phantasm ELF SHA-256:
`e760f9618bbf679811062b293f56cbb577fdcefba1b58661f94d2ea1f3758b72`.
The commit hook passes all 71 native tests plus the full-roster memory gate.
Targeted two-board captures also verified the repaired 4→5→6→7 and
9→10→11→12 sequences with zero spills.

## Caveats

- Counter scopes include live flywheel/DMA ISR time.
- Preset buckets own their following transition and therefore include the
  through-clear choreography.
- The final runs use normal choreography, not the accelerated diagnostic
  cycle.
- Page alignment is intentional performance state. Re-run the two-board full
  cycle after compiler, linker, or nearby flash-section changes.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShaderBall`, and
`HS_PROFILE_WINDOW=16`; `tools/profile_one.sh` performs locked build, flash,
capture, marker checks, and artifact attestation.
