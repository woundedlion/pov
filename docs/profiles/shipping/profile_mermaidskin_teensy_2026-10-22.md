# MermaidSkin on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/mermaidskin_ship.log`, captured 2026-10-22 20:03 America/Los_Angeles on COM4.
Replaces `profile_mermaidskin_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MermaidSkin 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh MermaidSkin profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:88552, data:154520, headers:8832   free for files:1779712` / `RAM1: variables:313024, code:16968, padding:15800   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/mermaidskin_ship_aa61e787859b_467b6d9838f7`. Profile ELF SHA-256: `467b6d9838f7d1f470b6a9c07965332840f98ed8c648d69aa71d3819df09537e`.

## Frame cadence

**Pass aggregate**: peak frame render is **61.132 ms** (frame 118), mean live-frame render is **59.755 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 53.001 ms/f across complete post-setup windows; its worst window is 53.421 ms/f (frames 33–64).

Setup frame 1 rendered 103.699 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:59) recorded peak 🟢 39.18 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 97–128)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.391 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.137 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 53.309 | 1.00 | 85.44% |
| ↳ ↳ `fx_prepare_column` | 1.333 | 144.00 | 2.14% |
| ↳ `fx_prepare_frame` | 4.090 | 1.00 | 6.56% |
| ↳ `fx_advance` | 1.933 | 1.00 | 3.10% |
| ↳ `fx_timeline_step` | 0.125 | 1.00 | 0.20% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 2.342 | 1.00 | 3.75% |

Wall min/avg/max = 60.658/62.402/64.214 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2000454 us
isr_wake       n=36874 cyc min/avg/max=331/1225/17686 ns min/avg/max=551/2041/29476 total=75344 us cpu=3.76%
isr_pack       n=4609 cyc min/avg/max=5321/5770/7337 ns min/avg/max=8868/9616/12228 total=44327 us cpu=2.21%
isr_dma_submit n=4609 cyc min/avg/max=387/568/4739 ns min/avg/max=645/946/7898 total=4368 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 53.309 ms/f, 85.44% of the peak window.
2. `fx_prepare_frame` — 4.090 ms/f, 6.56% of the peak window.
3. `fx_advance` — 1.933 ms/f, 3.10% of the peak window.
4. `pov_preserve_half` — 0.137 ms/f, 0.22% of the peak window.
5. `fx_timeline_step` — 0.125 ms/f, 0.20% of the peak window.
6. `canvas_clear` — 0.098 ms/f, 0.16% of the peak window.

README cells: peak 🟢 61.132, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=MermaidSkin  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (MermaidSkin 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (184.3 ppm, frames 1-32)

INVALID: MermaidSkin
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (184.3 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MermaidSkin`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
