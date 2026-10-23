# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/displacementfield_ship.log`, captured 2026-10-22 19:02 America/Los_Angeles on COM4.
Replaces `profile_displacementfield_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 121 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh DisplacementField profile 121 32` |

Image size (`profile`, this effect only): `FLASH: code:117952, data:154896, headers:8752   free for files:1750016` / `RAM1: variables:313472, code:25544, padding:7224   free for local variables:178048` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/displacementfield_ship_aa61e787859b_d026cb280789`. Profile ELF SHA-256: `d026cb2807897257526e8e6cbacc550fab6d2c0a0cd2b388a5286ee09fc522c7`.

## Frame cadence

**Pass aggregate**: peak frame render is **51.350 ms** (frame 603), mean live-frame render is **41.606 ms**, and **0/1887** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`df_timeline_step` averages 41.407 ms/f across complete post-setup windows; its worst window is 49.084 ms/f (frames 321–352).

Setup frame 1 rendered 18.995 ms. Published live-frame figures use complete counter windows through frame 1888; a trailing partial window is excluded.

The previous shipping report (2026-10-06 19:24) recorded peak 🟢 33.66 and spilled 🟢 0/2527 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 577–608)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `ring_stack_table` MIXED-PARENT | 1.818 | 2.00 | 2.91% |
| `frame` | 62.476 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.146 | 1.00 | 0.23% |
| ↳ `df_timeline_step` | 48.368 | 1.00 | 77.42% |
| ↳ ↳ `df_draw_rings` | 48.169 | 1.00 | 77.10% |
| ↳ ↳ ↳ `df_fused_scan` | 25.021 | 1.00 | 40.05% |
| ↳ ↳ ↳ `df_lut_bake` | 16.551 | 42.59 | 26.49% |
| ↳ ↳ ↳ ↳ `df_hue_table_prep` | 1.078 | 21.16 | 1.73% |
| ↳ ↳ ↳ ↳ `df_octave_noise` | 11.362 | 42.59 | 18.19% |
| ↳ ↳ ↳ `df_chunk_cull` | 1.573 | 43.91 | 2.52% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 13.823 | 1.00 | 22.13% |
| ↳ `df_prepare_fields` | 0.002 | 1.00 | 0.00% |

Wall min/avg/max = 60.535/62.479/64.624 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `df_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 1 | — | 51.350 | 0/898 | 49.084 | 27/28 |
| 2 | — | 43.390 | 0/989 | 41.284 | 30/30 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2002990 us
isr_wake       n=36921 cyc min/avg/max=364/1169/18523 ns min/avg/max=606/1948/30871 total=71962 us cpu=3.59%
isr_pack       n=4615 cyc min/avg/max=4475/5267/6987 ns min/avg/max=7458/8778/11645 total=40517 us cpu=2.02%
isr_dma_submit n=4615 cyc min/avg/max=363/564/4159 ns min/avg/max=605/940/6931 total=4343 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `df_timeline_step` — 48.368 ms/f, 77.42% of the peak window.
2. `pov_preserve_half` — 0.146 ms/f, 0.23% of the peak window.
3. `canvas_clear` — 0.095 ms/f, 0.15% of the peak window.
4. `df_prepare_fields` — 0.002 ms/f, 0.00% of the peak window.

README cells: peak 🟢 51.350 (2), spilled 🟢 0/1887 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=DisplacementField  windows=59
  [PASS] captured >=3 windows (59)
  [PASS] every window names one effect and resolution (DisplacementField 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] cycle returns to its first index (peak idx 2, 2 distinct)
  [PASS] all 2 presets visited (2)
  [PASS] per-frame render is render, not wall (0 of 59 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 59 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): df_chunk_cull, df_fused_scan, df_hue_table_prep, df_lut_bake, ring_stack_table
  [PASS] scope 'df_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (46.0 ppm, frames 321-352)

INVALID: DisplacementField
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (46.0 ppm, frames 321-352). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
