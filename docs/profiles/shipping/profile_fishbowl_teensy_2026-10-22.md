# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/fishbowl_ship.log`, captured 2026-10-22 18:59 America/Los_Angeles on COM4.
Replaces `profile_fishbowl_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Fishbowl 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Fishbowl profile 60 32` |

Image size (`profile`, this effect only): `FLASH: code:112032, data:152928, headers:8448   free for files:1758208` / `RAM1: variables:313024, code:34008, padding:31528   free for local variables:145728` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/fishbowl_ship_aa61e787859b_9545adc18572`. Profile ELF SHA-256: `9545adc1857259f770097540d19f4e403d0440c788e33ab3fc03657e722ca952`.

## Frame cadence

**Pass aggregate**: peak frame render is **28.480 ms** (frame 215), mean live-frame render is **22.762 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fish_build_vertices` averages 13.113 ms/f across complete post-setup windows; its worst window is 13.583 ms/f (frames 481–512).

Setup frame 1 rendered 1.779 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

The previous shipping report (2026-10-06 19:07) recorded peak 🟢 23.78 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 193–224)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.398 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.142 | 1.00 | 0.23% |
| ↳ `fish_multiline_draw` | 10.982 | 1.00 | 17.60% |
| ↳ `fish_build_vertices` | 13.526 | 1.00 | 21.68% |
| ↳ `fish_noise_prepare` | 0.002 | 1.00 | 0.00% |
| ↳ `fish_timeline_step` | 0.239 | 1.00 | 0.38% |
| ↳ ↳ `animation_noise_clock` | 0.000 | 1.00 | 0.00% |
| ↳ `canvas_clear` | 0.090 | 1.00 | 0.14% |
| ↳ `canvas_buffer_wait` | 36.742 | 1.00 | 58.88% |

Wall min/avg/max = 57.893/62.406/66.763 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2000848 us
isr_wake       n=36882 cyc min/avg/max=333/1499/21589 ns min/avg/max=555/2498/35981 total=92159 us cpu=4.60%
isr_pack       n=4610 cyc min/avg/max=3614/4047/5601 ns min/avg/max=6023/6745/9335 total=31094 us cpu=1.55%
isr_dma_submit n=4610 cyc min/avg/max=420/561/1547 ns min/avg/max=700/935/2578 total=4314 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fish_build_vertices` — 13.526 ms/f, 21.68% of the peak window.
2. `fish_multiline_draw` — 10.982 ms/f, 17.60% of the peak window.
3. `fish_timeline_step` — 0.239 ms/f, 0.38% of the peak window.
4. `pov_preserve_half` — 0.142 ms/f, 0.23% of the peak window.
5. `canvas_clear` — 0.090 ms/f, 0.14% of the peak window.
6. `fish_noise_prepare` — 0.002 ms/f, 0.00% of the peak window.

README cells: peak 🟢 28.480, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=Fishbowl  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (Fishbowl 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fish_build_vertices' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (107.4 ppm, frames 481-512)

INVALID: Fishbowl
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (107.4 ppm, frames 481-512). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Fishbowl`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
