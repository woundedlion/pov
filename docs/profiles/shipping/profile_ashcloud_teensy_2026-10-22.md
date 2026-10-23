# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/ashcloud_ship.log`, captured 2026-10-22 20:05 America/Los_Angeles on COM4.
Replaces `profile_ashcloud_teensy_2026-10-17.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh AshCloud profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:88504, data:154672, headers:8728   free for files:1779712` / `RAM1: variables:313024, code:16968, padding:15800   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/ashcloud_ship_aa61e787859b_7fdaa1be9eb8`. Profile ELF SHA-256: `7fdaa1be9eb83afcee92973076b9c0583e07d589f265b694591dcbd33518fbd0`.

## Frame cadence

**Pass aggregate**: peak frame render is **46.831 ms** (frame 864), mean live-frame render is **44.116 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 37.810 ms/f across complete post-setup windows; its worst window is 39.771 ms/f (frames 545–576).

Setup frame 1 rendered 71.679 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-17 17:12) recorded peak 🟢 46.27 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 833–864)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.410 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 38.964 | 1.00 | 62.43% |
| ↳ ↳ `fx_prepare_column` | 1.148 | 144.00 | 1.84% |
| ↳ `fx_prepare_frame` | 3.691 | 1.00 | 5.91% |
| ↳ `fx_advance` | 1.948 | 1.00 | 3.12% |
| ↳ `fx_timeline_step` | 0.115 | 1.00 | 0.18% |
| ↳ `canvas_clear` | 0.091 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 17.109 | 1.00 | 27.41% |

Wall min/avg/max = 59.965/62.423/64.730 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2000821 us
isr_wake       n=36881 cyc min/avg/max=328/1226/14523 ns min/avg/max=546/2043/24205 total=75419 us cpu=3.76%
isr_pack       n=4610 cyc min/avg/max=5328/5718/7344 ns min/avg/max=8880/9530/12240 total=43939 us cpu=2.19%
isr_dma_submit n=4610 cyc min/avg/max=391/564/4750 ns min/avg/max=651/940/7916 total=4337 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 38.964 ms/f, 62.43% of the peak window.
2. `fx_prepare_frame` — 3.691 ms/f, 5.91% of the peak window.
3. `fx_advance` — 1.948 ms/f, 3.12% of the peak window.
4. `pov_preserve_half` — 0.141 ms/f, 0.23% of the peak window.
5. `fx_timeline_step` — 0.115 ms/f, 0.18% of the peak window.
6. `canvas_clear` — 0.091 ms/f, 0.15% of the peak window.

README cells: peak 🟢 46.831, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

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
  [FAIL] root cyc/600 vs wall sum within 5 ppm (214.4 ppm, frames 545-576)

INVALID: AshCloud
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (214.4 ppm, frames 545-576). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
