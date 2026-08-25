# DisplacementField on-device profile - Teensy 4.0, segmented mode (2026-08-18, **selective -O3**)

Point-in-time snapshot (regenerate with `bash tools/profile_one.sh DisplacementField profile 90 32`).
Raw capture: `build/prof/displacementfield_ship.log`. Replaces
`profile_displacementfield_teensy_2026-07-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: `-Os` + newlib-nano + DMA LEDs + `HS_PROFILE_ENABLE`; distorted ring-stack scan, SDF distance chain, and `draw_rings` use selective `HS_O3` |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288x144, single-entry playlist, tip `d0514593e` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile 70 32` |

Image size: `FLASH: code:75632, data:149548, headers:8292` / `RAM1:
variables:315232, code:40360, padding:25176, free:143520` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 449-480 root counter cycles / 600 MHz
matches the measured wall sum within **1.1 ppm**
(`1,500,494,018 / 600 = 2,500,823.36 us` vs wall `2,500,826 us`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer):
`df_timeline_step` avg 47.66 ms/f, peak frame render **59.84 ms**
(frames 577-608), spilled **0/1408 frames (0%)**.

The detailed stage hierarchy below is the pre-fix regression capture, retained
as attribution evidence. Commit `6f73b10ec` moved the fused per-column ring
scan lambda from ITCM to cached flash; `5f985b55f` restores it to ITCM. The
larger remaining regression came from `b3b5c5bc5`: shrinking ColorWipe removed
an unused initial target-palette construction and shifted the shared RNG stream
before noise, orientation, and ball seeding. Commit `d0514593e` reserves that
draw without restoring the discarded palette object. The full Phantasm build
uses 188,680 B of RAM1 code, leaving 4,856 B against the hard ceiling and
passing the 3,072 B headroom ratchet.

Temporary deep scopes identify the remaining peak work, but perturb cadence:
`df_scan_distance` about 32.67 ms/f, LUT field bake 8.60 ms/f, shared scan
frame setup 5.18 ms/f, plot 3.25 ms/f, LUT hue bake 2.35 ms/f, hue-table
rotation 2.00 ms/f, ring shader 1.82 ms/f, and shape build 1.02 ms/f.

A display window is 62.5 ms; the effect renders one quadrant, about 10,368 px.
The pass mostly holds 16 fps; 16 frames in two adjacent peak windows exceed
62.5 ms. The `canvas_buffer_wait` scope is the round-up idle to the next
display flip and is excluded from render.

## Phase-by-phase readout

Phase schedule: one continuous regime where ring population rises and falls
through the breath cycle.

### Dwell peak - worst of capture (frames 481-512)

```text
frame                  78.03 ms  46.82 Mcyc 100%
  df_timeline_step     61.00 ms  36.60 Mcyc  78%
    df_draw_rings      60.91 ms  36.54 Mcyc  78%
      df_hue_table_prep 1.81 ms   1.08 Mcyc   2%  x22.2  81 us/call
      df_lut_bake      10.91 ms   6.54 Mcyc  14%  x45.2 241 us/call
      df_chunk_cull     1.08 ms 648.60 kcyc   1%  x46.9  23 us/call
      df_fused_scan    45.65 ms  27.39 Mcyc  58%  x1.0
        filter_blend    1.09 ms 651.42 kcyc   1%  x10414 63 cyc/blend
  canvas_clear         89.69 us  53.83 kcyc   0%
  canvas_buffer_wait   16.94 ms  10.17 Mcyc  22%
  df_prepare_fields     0.06 us      56 cyc   0%
```

Wall min/avg/max = 59.442/78.031/125.302 ms. `df_fused_scan` dominates
this peak regime at **45.65 ms/f**, while `df_lut_bake` adds 10.91 ms/f.

### Dwell trough (frames 897-928)

```text
frame                  63.00 ms  37.80 Mcyc 100%
  df_timeline_step     21.25 ms  12.75 Mcyc  34%
    df_draw_rings      21.14 ms  12.69 Mcyc  34%
      df_hue_table_prep 0.44 ms 261.14 kcyc   1%  x10.5  41 us/call
      df_lut_bake       1.96 ms   1.18 Mcyc   3%  x12.2 160 us/call
      df_chunk_cull     0.32 ms 194.44 kcyc   1%  x13.1  25 us/call
      df_fused_scan    17.41 ms  10.45 Mcyc  28%  x1.0
        filter_blend    0.92 ms 551.79 kcyc   1%  x9157  60 cyc/blend
  canvas_clear         86.00 us  51.62 kcyc   0%
  canvas_buffer_wait   41.65 ms  24.99 Mcyc  66%
  df_prepare_fields     1.28 us     783 cyc   0%
```

Wall min/avg/max = 58.451/62.996/66.661 ms. `df_fused_scan` drops to
17.41 ms/f and `df_lut_bake` to 1.96 ms/f, so this is the same regime
under lower dwell pressure.

### Per-pixel figures

Peak: 10,414 blended px/frame against 10,368 px = **1.00x** coverage;
`filter_blend` costs 62.6 cyc/blend and `df_fused_scan` consumes about
2,629.7 scan cycles per blended pixel.

Trough: 9,157 blended px/frame = **0.88x** coverage; `filter_blend` costs
60.3 cyc/blend and `df_fused_scan` consumes about 1,141.9 scan cycles per
blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1439/frame  min/avg/max 0.705/1.986/12.176 us  cpu 3.66%
isr_pack         180/frame  min/avg/max 6.368/7.013/9.451 us   cpu 1.61%
isr_dma_submit   180/frame  min/avg/max 0.608/0.940/2.933 us   cpu 0.21%
```

- Submit costs about 13% as much CPU as pack; wire transfer is asynchronous.
- ISR share is about 5.48%, leaving about 59.1 ms of CPU budget in a 62.5 ms
  window. The 59.84 ms peak needs about a 1.3% reduction against that budget.

## Summary ranking

1. `df_fused_scan` - 58% of the wall frame, 45.65 ms.
2. `df_lut_bake` - 14%, 10.91 ms.
3. `df_hue_table_prep` - 2%, 1.81 ms.
4. `df_chunk_cull` - 1%, 1.08 ms.

No matching WASM/native DisplacementField entry was available in the perf
ledger for a direct comparison.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under `df_fused_scan`; its subtree is hidden in
  windows where that parent has no calls.
- Per-pixel profiling is reported from the existing filter counter; no extra
  per-pixel scope was added for this capture.
- This is the shipping selective-`O3` image; distorted ring-stack scan, the
  SDF distance chain, and `draw_rings` cross `HS_O3` regions.
- No dwell-compression knobs were used.
- Captured from tip `74851e0e8` with the palette callback that landed as
  `bde19bb72`. A clean-parent recapture reproduced 63.06 ms and 16/1056
  spills, so the spill predates the callback.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`,
`HS_PROFILE_WINDOW=32`; `just profile DisplacementField` builds, flashes,
and captures under the device lock.
