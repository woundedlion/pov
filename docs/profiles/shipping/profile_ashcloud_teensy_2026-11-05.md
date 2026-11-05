# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AshCloud`).
Raw capture: `C:/work/temp/Holosphere-20261105-ashcloud-fix-notes/fixed.log`. Captured 2026-11-05 08:38 America/Los_Angeles on COM3. Replaces the 00:10 capture in this report; its original text is retained at `C:/work/temp/Holosphere-20261105-ashcloud-fix-notes/previous_good_profile.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3 serial 6362410, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths; shared AshCloud shader in ITCM |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AshCloud, 288×144, single-entry playlist, measured source `b4ae7a949f9819698f0a9c10a32df5c6df683183` |
| Method | 60 s, 32-frame counter windows; default one-hour epoch, no dwell compression |
| Reproduce | `HS_PROFILE_TREE=C:/work/Holosphere HS_TEENSY_PORT=COM3 bash tools/profile_one.sh AshCloud profile 60 32` |

Profile image size:

```text
teensy_size:   FLASH: code:121432, data:182688, headers:8200   free for files:1719296
teensy_size:    RAM1: variables:314144, code:40072, padding:25464   free for local variables:144608
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Shipping Phantasm image size:

```text
teensy_size:   FLASH: code:1177320, data:791364, headers:8660   free for files:54272
teensy_size:    RAM1: variables:313760, code:192248, padding:4360   free for local variables:13920
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Against the matched regressed source `8e86fcd1be02056fcd5f7637b5a1cc038fa728dc`, Phantasm RAM1 code is 192,248 B (+96), RAM1 variables 313,760 B (+0), and FLASH data 791,364 B (+0). The shipping memory/layout gate passes with 13,920 B stack headroom and 4,360 B remaining within six ITCM banks. AshCloud sources landed unchanged as `b3391d0f4` after unrelated mesh commits were integrated.

Profile ELF SHA-256: `1bd95022821e39160607000eda23fb7a6d9fd46cd21b94361e0e44f23bca28f9`. Phantasm ELF SHA-256: `b9aa6e131df11baa07bd765aaabb13d301ebc79e63ff22564eb5abeda1fb39af`. Both archived ELFs were rehashed while generating this report; adjacent provenance and envdump files retain compiler/ABI and flags.

Exactness cross-check: frames 897–928, root cyc/600 versus measured wall sum differs by **0.774 ppm**, within 5 ppm. Parser validation passes.

## Frame cadence

**Steady render, frame 1 excluded:** peak **58.057 ms** (frame 920), **0/927 spilled** (0.00%). Setup frame 1: 80.859 ms. Complete windows cover frames 1–928; trailing incomplete rows are excluded.

Steady mean render is 50.763 ms/frame. `fx_shader_draw` averages 47.586 ms/frame across complete windows; its worst window is 50.796 ms/frame (frames 897–928). A display window is 62.5 ms at 480 RPM; one quadrant is 10,368 pixels. `canvas_buffer_wait` is the round-up idle to the next display flip. All measured steady frames fit 16 fps, with 4.443 ms peak margin.

Matched before/after captures use the same COM3 board and stock harness. Both raw logs pass validation. Setup frame 1 is excluded below.

| Capture | Frames | Peak ms | Spilled | Mean render ms | Shader mean ms |
|---|--:|--:|--:|--:|--:|
| Regressed `8e86fcd1b` | 703 | 66.104 | 102/703 (14.51%) | 59.259 | 56.079 |
| Fixed, matched frames 2–704 | 703 | 57.282 | 0/703 (0.00%) | 50.114 | — |
| Fixed, complete capture | 927 | 58.057 | 0/927 (0.00%) | 50.763 | 47.586 |

Regressed log: `C:/work/temp/Holosphere-20261105-ashcloud-fix-notes/baseline.log`. Matched mean render falls 15.43%. The earlier 00:10 COM4 capture at `e00f48da5` recorded 56.969 ms peak and 0/927 spills; it used factored camera arithmetic and is historical context rather than the matched control.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. Scopes are inclusive; nested scopes overlap and are not additive.

### Startup (frames 1–32)

```text
frame                                 62.389 ms 37.434 Mcyc 100.0%
  pov_preserve_half                    0.136 ms  0.082 Mcyc   0.2%
    x0.97/frame; 140.581 us/call
  fx_shader_draw                      42.980 ms 25.788 Mcyc  68.9%
    x1.00/frame; 42980.250 us/call
    fx_prepare_column                  3.153 ms  1.892 Mcyc   5.1%
      x144.00/frame; 21.893 us/call
  fx_palette_base_bake                 0.174 ms  0.104 Mcyc   0.3%
    x1.00/frame; 173.688 us/call
    fx_palette_rgb_bake                0.054 ms  0.032 Mcyc   0.1%
      x1.00/frame; 53.781 us/call
  fx_prepare_frame                     0.917 ms  0.550 Mcyc   1.5%
    x1.00/frame; 917.000 us/call
  fx_advance                           1.119 ms  0.671 Mcyc   1.8%
    x1.00/frame; 1118.938 us/call
    palette_cycle_step                 1.019 ms  0.611 Mcyc   1.6%
      x1.00/frame; 1018.812 us/call
  fx_timeline_step                     0.121 ms  0.072 Mcyc   0.2%
    x1.00/frame; 120.625 us/call
  canvas_clear                         0.100 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.062 us/call
  canvas_buffer_wait                  16.571 ms  9.943 Mcyc  26.6%
    x1.00/frame; 16571.438 us/call
```

Wall min/avg/max: 39.319/62.389/80.862 ms. Window mean render is 45.818 ms. The startup frame includes initialization; display synchronization contributes to wall time.

### Peak-render regime (frames 897–928)

```text
frame                                 62.444 ms 37.467 Mcyc 100.0%
  pov_preserve_half                    0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 137.000 us/call
  fx_shader_draw                      50.796 ms 30.477 Mcyc  81.3%
    x1.00/frame; 50795.781 us/call
    fx_prepare_column                  3.494 ms  2.096 Mcyc   5.6%
      x144.00/frame; 24.262 us/call
  fx_palette_base_bake                 0.614 ms  0.369 Mcyc   1.0%
    x1.00/frame; 614.406 us/call
    fx_palette_rgb_bake                0.486 ms  0.291 Mcyc   0.8%
      x1.00/frame; 485.688 us/call
  fx_prepare_frame                     1.091 ms  0.654 Mcyc   1.7%
    x1.00/frame; 1090.812 us/call
  fx_advance                           1.149 ms  0.689 Mcyc   1.8%
    x1.00/frame; 1148.875 us/call
    palette_cycle_step                 1.093 ms  0.656 Mcyc   1.8%
      x1.00/frame; 1092.906 us/call
  fx_timeline_step                     0.128 ms  0.077 Mcyc   0.2%
    x1.00/frame; 128.375 us/call
  canvas_clear                         0.097 ms  0.058 Mcyc   0.2%
    x1.00/frame; 97.344 us/call
  canvas_buffer_wait                   8.169 ms  4.902 Mcyc  13.1%
    x1.00/frame; 8169.469 us/call
```

Wall min/avg/max: 57.522/62.444/67.359 ms. Window mean render is 54.275 ms. This window contains the steady peak and the largest shader mean; render remains below the display window.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                 62.398 ms 37.439 Mcyc 100.0%
  pov_preserve_half                    0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 140.156 us/call
  fx_shader_draw                      44.312 ms 26.587 Mcyc  71.0%
    x1.00/frame; 44312.062 us/call
    fx_prepare_column                  2.724 ms  1.634 Mcyc   4.4%
      x144.00/frame; 18.913 us/call
  fx_palette_base_bake                 0.200 ms  0.120 Mcyc   0.3%
    x1.00/frame; 199.750 us/call
    fx_palette_rgb_bake                0.067 ms  0.040 Mcyc   0.1%
      x1.00/frame; 66.844 us/call
  fx_prepare_frame                     1.012 ms  0.607 Mcyc   1.6%
    x1.00/frame; 1011.531 us/call
  fx_advance                           1.139 ms  0.684 Mcyc   1.8%
    x1.00/frame; 1139.219 us/call
    palette_cycle_step                 1.096 ms  0.658 Mcyc   1.8%
      x1.00/frame; 1096.000 us/call
  fx_timeline_step                     0.123 ms  0.074 Mcyc   0.2%
    x1.00/frame; 123.156 us/call
  canvas_clear                         0.095 ms  0.057 Mcyc   0.2%
    x1.00/frame; 95.094 us/call
  canvas_buffer_wait                  15.117 ms  9.070 Mcyc  24.2%
    x1.00/frame; 15117.188 us/call
```

Wall min/avg/max: 58.554/62.398/66.192 ms. Window mean render is 47.281 ms. Lower shader work leaves more synchronization idle; the sampling resolution and cadence are unchanged.

### Per-pixel figures

Per-pixel blend counters are disabled in this shipping timing image. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake            1153.84/frame
  min/avg/max 0.535/1.958/26.555 us; cpu 3.61%
  isr_pack           144.22/frame
  min/avg/max 7.551/8.988/12.170 us; cpu 2.07%
  isr_dma_submit     144.22/frame
  min/avg/max 0.668/0.915/6.938 us; cpu 0.21%
```

Raw ISR readout and separately measured interval:

```text
ISR window=2003060 us
isr_wake       n=36923 cyc min/avg/max=321/1175/15933 ns min/avg/max=535/1958/26555 total=72333 us cpu=3.61%
isr_pack       n=4615 cyc min/avg/max=4531/5393/7302 ns min/avg/max=7551/8988/12170 total=41485 us cpu=2.07%
isr_dma_submit n=4615 cyc min/avg/max=401/549/4163 ns min/avg/max=668/915/6938 total=4226 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. The inclusive wake share of 3.61% leaves approximately 60.244 ms of foreground CPU per 62.5 ms window. Render scopes already include ISR time, so no second ISR subtraction is needed; measured steady rendering needs no speedup to fit one window.

## Summary ranking

1. `fx_shader_draw` — 50.796 ms/frame, 81.35% of this root window.
2. `fx_prepare_column` — 3.494 ms/frame, 5.60% of this root window.
3. `fx_advance` — 1.149 ms/frame, 1.84% of this root window.
4. `palette_cycle_step` — 1.093 ms/frame, 1.75% of this root window.
5. `fx_prepare_frame` — 1.091 ms/frame, 1.75% of this root window.

Commit `04e7d2166` disabled unqualified camera factoring but also made `!factored` bypass the shared color shader. The earlier PC profile (`C:/work/temp/ashcloud-20261030/ashcloud-final-pc_report.txt`) attributed 57.7% of render samples to the ITCM shared shader. The fix admits shared color independently and executes canonical per-pixel rotation inside that same ITCM kernel. Native/WASM execution is a correctness check, not a hardware timing substitute.

Debug qualification: 42,614 assertions passed, zero failed. Optimized `-O3 -DNDEBUG -ffast-math -fno-finite-math-only` qualification: 46,453 assertions passed, zero failed. Both inspect 1,380,096 samples. Strict output is bitwise equal; optimized maxima are one RGB16 code, one premultiplied code, zero alpha error and zero camera-direction error. Replacing selected shared output with transparent black fails the column regression assertions (7,680 failures).

README cells: peak 🟢 58.057, spilled 🟢 0/927 (0.00%).

README setup frame 1: 80.859 ms. Captured 2026-11-05 08:38 America/Los_Angeles.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled; per-pixel counters are disabled here.
- Shipping selective `HS_O3` shader/draw paths are active; no dwell-compression knobs were used.
- Capture builds began with the isolated source changes uncommitted; the exact source was committed as `b4ae7a949f9819698f0a9c10a32df5c6df683183` before ELF attestation and capture. Archived binaries and hashes identify the measured images.
- The 60 s pass bounds the measured default-control interval; it is not an exhaustive parameter sweep.

Parser validation:

```text
effect=AshCloud  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (AshCloud 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.8 ppm, frames 897-928)

VALID: AshCloud
```

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`. `just profile AshCloud 60` builds, flashes and captures through the shared-device lock.
