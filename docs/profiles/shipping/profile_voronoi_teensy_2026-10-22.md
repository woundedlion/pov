# Voronoi on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/voronoi_ship.log`, captured 2026-10-22 19:13 America/Los_Angeles on COM4.
Replaces `profile_voronoi_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Voronoi 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 32-frame windows, 30 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Voronoi profile 30 32` |

Image size (`profile`, this effect only): `FLASH: code:50216, data:147068, headers:8540   free for files:1825792` / `RAM1: variables:312960, code:10344, padding:22424   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/voronoi_ship_aa61e787859b_1d4150382136`. Profile ELF SHA-256: `1d4150382136ab034e747535816627d773f294e8549f12d3032d7ee82eed6ab7`.

## Frame cadence

**Pass aggregate**: peak frame render is **9.243 ms** (frame 2), mean live-frame render is **8.425 ms**, and **0/447** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`vo_shade` averages 7.758 ms/f across complete post-setup windows; its worst window is 7.855 ms/f (frames 417–448).

Setup frame 1 rendered 13.235 ms. Published live-frame figures use complete counter windows through frame 448; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:25) recorded peak 🟢 8.00 and spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 1–32)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 59.237 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.148 | 0.97 | 0.25% |
| ↳ `vo_shade` DUPLICATE-NAME | 0.393 | 0.03 | 0.66% |
| ↳ `vo_kdtree` | 0.344 | 1.00 | 0.58% |
| ↳ `vo_animate` | 0.043 | 1.00 | 0.07% |
| ↳ `canvas_clear` | 0.089 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 50.597 | 1.00 | 85.41% |

Wall min/avg/max = 9.244/59.251/62.931 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=2382796 us
isr_wake       n=43922 cyc min/avg/max=268/990/8578 ns min/avg/max=446/1650/14296 total=72517 us cpu=3.04%
isr_pack       n=4341 cyc min/avg/max=5244/5531/7179 ns min/avg/max=8740/9218/11965 total=40020 us cpu=1.67%
isr_dma_submit n=4341 cyc min/avg/max=375/561/620 ns min/avg/max=625/935/1033 total=4062 us cpu=0.17%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `vo_shade` — 0.393 ms/f, 0.66% of the peak window.
2. `vo_kdtree` — 0.344 ms/f, 0.58% of the peak window.
3. `pov_preserve_half` — 0.148 ms/f, 0.25% of the peak window.
4. `canvas_clear` — 0.089 ms/f, 0.15% of the peak window.
5. `vo_animate` — 0.043 ms/f, 0.07% of the peak window.

README cells: peak 🟢 9.243, spilled 🟢 0/447 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=Voronoi  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (Voronoi 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): vo_shade
  [PASS] scope 'vo_shade' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (244.3 ppm, frames 417-448)

INVALID: Voronoi
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (244.3 ppm, frames 417-448). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Voronoi`, `HS_PROFILE_WINDOW=32`; the Reproduce command builds, flashes and captures under the COM4 device lock.
