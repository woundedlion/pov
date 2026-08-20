# DisplacementField on-device profile - Teensy 4.0, segmented mode (2026-08-18, **-O3**)

Global-`O3` twin of the
[shipping capture](../shipping/profile_displacementfield_teensy_2026-08-18.md).
Raw capture: `build/prof/displacementfield_o3.log`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect reference ceiling |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288x144, single-entry playlist, tip `d0514593e` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile_o3 90 32` |

Image size: `FLASH: code:100848, data:149124, headers:9100` / `RAM1:
variables:315264, code:61960, padding:3576, free:143488` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 449-480 root counter cycles / 600 MHz
matches the measured wall sum within **2.3 ppm**
(`1,200,804,795 / 600 = 2,001,341.33 us` vs wall `2,001,346 us`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer):
`df_timeline_step` avg 46.77 ms/f, peak frame render **58.96 ms**
(frames 577-608), spilled **0/1408 frames (0%)**.

The detailed stage hierarchy below is the pre-fix regression capture, retained
as attribution evidence. Restoring the fused per-column scan lambda to ITCM in
`5f985b55f` and restoring the geometry RNG stream in `d0514593e` bring global
`-O3` back below the ISR-adjusted frame budget with zero spills.

A display window is 62.5 ms; the effect renders one quadrant, about 10,368 px.
The pass holds 16 fps without a render spill. The `canvas_buffer_wait` scope
is display-sync idle and is excluded from render.

## Phase-by-phase readout

Phase schedule: one continuous regime matching shipping, with ring population
rising and falling through the breath cycle.

### Dwell peak - worst of capture (frames 481-512)

```text
frame                  62.43 ms  37.46 Mcyc 100%
  df_timeline_step     60.01 ms  36.00 Mcyc  96%
    df_draw_rings      59.92 ms  35.95 Mcyc  96%
      df_hue_table_prep 1.56 ms 937.01 kcyc   3%  x21.3  73 us/call
      df_lut_bake      10.30 ms   6.18 Mcyc  16%  x45.2 228 us/call
      df_chunk_cull     1.02 ms 612.43 kcyc   2%  x46.5  22 us/call
      df_fused_scan    45.63 ms  27.38 Mcyc  73%  x1.0
        filter_blend    1.12 ms 670.30 kcyc   2%  x10443 64 cyc/blend
  canvas_clear         90.22 us  54.14 kcyc   0%
  canvas_buffer_wait    2.33 ms   1.40 Mcyc   4%
  df_prepare_fields     0.09 us      68 cyc   0%
```

Wall min/avg/max = 58.707/62.425/65.953 ms. Global `-O3` leaves the fused
scan nearly unchanged at 45.63 ms/f; most of its gain over shipping is in the
supporting ring work.

### Dwell trough (frames 897-928)

```text
frame                  63.00 ms  37.80 Mcyc 100%
  df_timeline_step     22.41 ms  13.44 Mcyc  36%
    df_draw_rings      22.30 ms  13.38 Mcyc  35%
      df_hue_table_prep 0.47 ms 282.33 kcyc   1%  x13.8  34 us/call
      df_lut_bake       2.20 ms   1.32 Mcyc   3%  x17.6 125 us/call
      df_chunk_cull     0.42 ms 252.20 kcyc   1%  x17.6  24 us/call
      df_fused_scan    18.14 ms  10.88 Mcyc  29%  x1.0
        filter_blend    0.92 ms 550.99 kcyc   1%  x9152  60 cyc/blend
  canvas_clear         85.88 us  51.53 kcyc   0%
  canvas_buffer_wait   40.50 ms  24.30 Mcyc  64%
  df_prepare_fields     1.16 us     702 cyc   0%
```

Wall min/avg/max = 58.462/63.000/65.525 ms. The trough remains the same
structural regime, with lower ring count and dwell pressure.

### Per-pixel figures

Peak: 10,443 blended px/frame against 10,368 px = **1.01x** coverage;
`filter_blend` costs 64.2 cyc/blend and `df_fused_scan` consumes about
2,622.0 scan cycles per blended pixel.

Trough: 9,152 blended px/frame = **0.88x** coverage; `filter_blend` costs
60.2 cyc/blend and `df_fused_scan` consumes about 1,189.0 scan cycles per
blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1151/frame  min/avg/max 0.575/1.870/15.691 us  cpu 3.44%
isr_pack         144/frame  min/avg/max 6.198/6.921/13.863 us  cpu 1.59%
isr_dma_submit   144/frame  min/avg/max 0.673/0.936/2.206 us   cpu 0.21%
```

- Submit costs about 14% as much CPU as pack; wire transfer is asynchronous.
- ISR share is about 5.24%, leaving about 59.2 ms of CPU budget in a 62.5 ms
  window. The 58.96 ms peak is inside that budget.

## Summary ranking

1. `df_fused_scan` - 73% of the wall frame, 45.63 ms.
2. `df_lut_bake` - 16%, 10.30 ms.
3. `df_hue_table_prep` - 3%, 1.56 ms.
4. `df_chunk_cull` - 2%, 1.02 ms.

No matching WASM/native DisplacementField entry was available in the perf
ledger for a direct comparison.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under `df_fused_scan`; its subtree is hidden in
  windows where that parent has no calls.
- Per-pixel profiling is reported from the existing filter counter; no extra
  per-pixel scope was added for this capture.
- Global `-O3` is a single-effect ceiling and is not shippable at roster
  scale.
- No dwell-compression knobs were used.
- Captured from tip `74851e0e8` with the palette callback that landed as
  `bde19bb72`.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`,
`HS_PROFILE_WINDOW=32`; use
`bash tools/profile_one.sh DisplacementField profile_o3 90 32`.

## Global `-O3` vs selective `-O3`

Global `-O3` lowers the peak from **59.84 ms to 58.96 ms**
(**0.66 ms, 1.0%**) and adds **25,216 B FLASH** and **21,600 B RAM1 code**.
The shipping spill therefore comes from pre-existing workload headroom, not
the palette-mapping callback or a stale capture.
