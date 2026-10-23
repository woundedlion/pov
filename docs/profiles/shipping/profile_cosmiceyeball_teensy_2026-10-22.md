# CosmicEyeball on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/cosmiceyeball_ship.log`, captured 2026-10-22 20:21 America/Los_Angeles on COM4.
Replaces `profile_cosmiceyeball_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | CosmicEyeball 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh CosmicEyeball profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:85456, data:154328, headers:9048   free for files:1782784` / `RAM1: variables:313024, code:16536, padding:16232   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/cosmiceyeball_ship_aa61e787859b_b95964fc2fce`. Profile ELF SHA-256: `b95964fc2fcefaf47f59f14f468fa8cb821ccc4db2a4d263ffa50091c4d1ebe2`.

## Frame cadence

**Pass aggregate**: peak frame render is **28.921 ms** (frame 246), mean live-frame render is **28.075 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 24.856 ms/f across complete post-setup windows; its worst window is 24.912 ms/f (frames 513–544).

Setup frame 1 rendered 46.910 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-06 19:06) recorded peak 🟢 23.63 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 225–256)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.398 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.146 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 24.888 | 1.00 | 39.89% |
| ↳ ↳ `fx_prepare_column` | 1.862 | 144.00 | 2.98% |
| ↳ `fx_prepare_frame` | 0.613 | 1.00 | 0.98% |
| ↳ `fx_advance` | 2.054 | 1.00 | 3.29% |
| ↳ `fx_timeline_step` | 0.094 | 1.00 | 0.15% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 34.325 | 1.00 | 55.01% |

Wall min/avg/max = 61.357/62.410/63.459 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=1999891 us
isr_wake       n=36865 cyc min/avg/max=365/1206/16277 ns min/avg/max=608/2010/27128 total=74123 us cpu=3.70%
isr_pack       n=4608 cyc min/avg/max=5268/5745/7222 ns min/avg/max=8780/9575/12036 total=44121 us cpu=2.20%
isr_dma_submit n=4608 cyc min/avg/max=389/562/5354 ns min/avg/max=648/936/8923 total=4320 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 24.888 ms/f, 39.89% of the peak window.
2. `fx_advance` — 2.054 ms/f, 3.29% of the peak window.
3. `fx_prepare_frame` — 0.613 ms/f, 0.98% of the peak window.
4. `pov_preserve_half` — 0.146 ms/f, 0.23% of the peak window.
5. `fx_timeline_step` — 0.094 ms/f, 0.15% of the peak window.
6. `canvas_clear` — 0.089 ms/f, 0.14% of the peak window.

README cells: peak 🟢 28.921, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=CosmicEyeball  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (CosmicEyeball 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (163.7 ppm, frames 1-32)

INVALID: CosmicEyeball
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (163.7 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=CosmicEyeball`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
