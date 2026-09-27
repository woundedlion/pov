# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-26, **selective -O3**)

Point-in-time snapshot of **Preset 2, hypercube-flight**. Replaces the historical 2026-08-26 report; this fixed-preset capture is not a full-roster or full-cycle result.

[Raw capture](../evidence/hyperlattice_preset2_2026-09-26/ship.txt), [image provenance](../evidence/hyperlattice_preset2_2026-09-26/ship.provenance), and [build/size evidence](../evidence/hyperlattice_preset2_2026-09-26/ship_build.txt).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3, flywheel + DMA ISRs live |
| Image | `profile`; -Os base with shipping selective-O3 helpers; deep instrumentation off |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144; clean source `6b50b3de66891418e076e1f079826c6e252f0859` |
| Method | 70 s, 32-frame windows; `HS_PROFILE_PRESET=1` fixes UI Preset 2; motion and palette continue; no automatic preset transitions |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice profile 70 32 "-D HS_PROFILE_PRESET=1"` |

Image size:

```text
FLASH: code:63,440, data:148,524, headers:8,196, free:1,811,456
RAM1: variables:314,944, code:20,824, padding:11,944, free:176,576
RAM2: variables:520,064, free:4,224
```

Exactness cross-check: frames 417–448, root 1,203,125,622 cycles / 600 MHz versus wall sum 2,005,217 µs: **3.81 ppm**. The untouched capture was checked with `tools/parse_profile.py` in validate mode.

## Frame cadence

Runtime frames **2–1088**: mean render **44.200 ms**, peak **50.095 ms**, spilled **0/1087 (0.00%)**. Mean wall time is 62.440 ms. Startup frame 1 renders in **60.086 ms** and is excluded from all runtime statistics.

Scope and ISR summaries use complete windows **33–1088**, excluding the startup-containing window. Undumped trailing frames are outside the reported range.

The live quadrant has 10,368 samples. One display window is 62.5 ms at 480 RPM; peak render leaves 12.405 ms of headroom. `canvas_buffer_wait` is alignment idle, not render work.

## Phase-by-phase readout

Preset 2 remains selected throughout. Changing view orientation, ray intersections, and the palette morph vary cost within this held preset; there is no preset-cycle wrap to measure.

### Highest-cost window (frames 417–448)

```text
frame                       62.663 ms 37.598 Mcyc 100.0%
  pov_preserve_half          0.140 ms  0.084 Mcyc   0.2% x1.0
  hl_shader_draw            44.317 ms 26.590 Mcyc  70.7% x1.0
  hl_timeline_step           0.005 ms  0.003 Mcyc   0.0% x1.0
  canvas_clear               0.085 ms  0.051 Mcyc   0.1% x1.0
  canvas_buffer_wait        16.330 ms  9.798 Mcyc  26.1% x1.0
```

Wall min/avg/max: 57.826/62.663/67.407 ms. Render averages 46.333 ms; shader evaluation is the dominant cost. The wait aligns publication to the next display flip.

### Lowest-cost window (frames 289–320)

```text
frame                       62.303 ms 37.382 Mcyc 100.0%
  pov_preserve_half          0.139 ms  0.083 Mcyc   0.2% x1.0
  hl_shader_draw            36.658 ms 21.995 Mcyc  58.8% x1.0
  hl_timeline_step           0.005 ms  0.003 Mcyc   0.0% x1.0
  canvas_clear               0.084 ms  0.051 Mcyc   0.1% x1.0
  canvas_buffer_wait        23.649 ms 14.190 Mcyc  38.0% x1.0
```

Wall min/avg/max: 58.224/62.302/66.625 ms. Render averages 38.653 ms; shader evaluation is the dominant cost. The wait aligns publication to the next display flip.

### Per-pixel figures

The shader writes directly; no `filter_blend` population is available. The highest-cost window spends 2564.6 shader cycles per configured sample. Its shader call costs 44316.6 µs.

## Column-ISR / DMA marshaling cost

Rates are calls/frame; min/average/max are per-call microseconds, followed by CPU share.

```text
isr_wake       1152.1/f  0.61/1.65/18.80 us  3.04%
isr_pack        144.0/f  6.23/6.71/9.72 us  1.55%
isr_dma_submit  144.0/f  0.60/0.94/8.34 us  0.22%
```

Pack and submit measure CPU marshaling; LED wire transfer is asynchronous. Total ISR share is 4.80%, equivalent to 3.002 ms per display window and approximately 59.498 ms left for foreground work. Render counters already include ISR interruptions; do not subtract them twice. No speedup is required for the observed 16 fps cadence.

## Summary ranking

1. `hl_shader_draw`: 42.235 ms/frame, 67.6% of root cycles across complete live windows.
2. Palette preparation and other unscoped foreground work occupy the remainder of render time.
3. `canvas_buffer_wait` is available display-alignment idle.

No matched native/WASM timing comparison was collected.

## Caveats

- This measures only Preset 2, not transitions or the worst case over every preset.
- All scopes include ISR time because CYCCNT free-runs.
- Deep per-pixel scopes are disabled; their overhead is absent. `filter_blend` parenting does not affect this direct-write path.
- Shipping uses the selective-O3 shader scan and color/gamut helpers; O3 globally changes optimization and floating-point flags.
- Fixed selection pauses choreography, not spatial motion or palette cycling. No dwell compression or epoch stretch was used.
- Source was clean. A local Python 3.12 virtual environment supplied the repository-pinned PlatformIO because the host Python 3.14 SSL DLL was unavailable.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=32`, and `HS_PROFILE_PRESET=1`. Use the locked reproduce command above; plain `just profile HyperLattice` does not pin Preset 2.
