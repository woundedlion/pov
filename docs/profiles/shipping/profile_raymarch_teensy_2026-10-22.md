# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/raymarch_ship.log`, captured 2026-10-22 19:10 America/Los_Angeles on COM4.
Replaces `profile_raymarch_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Raymarch 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Raymarch profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:112296, data:192288, headers:8760   free for files:1718272` / `RAM1: variables:312896, code:25032, padding:7736   free for local variables:178624` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/raymarch_ship_aa61e787859b_ae10453108ee`. Profile ELF SHA-256: `ae10453108ee7fa779a59ab24f3e9a0569bd93592e16e7d7a65ce9c8e2c74f48`.

## Frame cadence

**Pass aggregate**: peak frame render is **57.674 ms** (frame 10), mean live-frame render is **52.162 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`rm_shader_draw` averages 48.564 ms/f across complete post-setup windows; its worst window is 49.650 ms/f (frames 321–352).

Setup frame 1 rendered 103.857 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:17) recorded peak 🟢 51.95 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 1–32)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 63.289 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.129 | 0.97 | 0.20% |
| ↳ `rm_shader_draw` | 50.560 | 1.00 | 79.89% |
| ↳ `rm_timeline_step` | 0.421 | 1.00 | 0.67% |
| ↳ `canvas_clear` | 0.092 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 9.177 | 1.00 | 14.50% |

Wall min/avg/max = 53.078/63.291/103.859 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2421891 us
isr_wake       n=44643 cyc min/avg/max=278/940/9461 ns min/avg/max=463/1566/15768 total=69996 us cpu=2.89%
isr_pack       n=4432 cyc min/avg/max=3594/4741/6609 ns min/avg/max=5990/7901/11015 total=35022 us cpu=1.44%
isr_dma_submit n=4432 cyc min/avg/max=371/562/622 ns min/avg/max=618/936/1036 total=4155 us cpu=0.17%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `rm_shader_draw` — 50.560 ms/f, 79.89% of the peak window.
2. `rm_timeline_step` — 0.421 ms/f, 0.67% of the peak window.
3. `pov_preserve_half` — 0.129 ms/f, 0.20% of the peak window.
4. `canvas_clear` — 0.092 ms/f, 0.15% of the peak window.

README cells: peak 🟢 57.674, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=Raymarch  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (Raymarch 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'rm_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (36.9 ppm, frames 1-32)

INVALID: Raymarch
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (36.9 ppm, frames 1-32). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Raymarch`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
