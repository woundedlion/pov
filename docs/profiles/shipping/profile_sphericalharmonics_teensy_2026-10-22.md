# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/sphericalharmonics_ship.log`, captured 2026-10-22 19:28 America/Los_Angeles on COM4.
Replaces `profile_sphericalharmonics_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 96 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh SphericalHarmonics profile 96 16 "-D HS_PROFILE_ORDERED_CYCLE"` |

Image size (`profile`, this effect only): `FLASH: code:43008, data:147376, headers:8272   free for files:1832960` / `RAM1: variables:312896, code:11864, padding:20904   free for local variables:178624` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/sphericalharmonics_ship_aa61e787859b_a4a28d281491`. Profile ELF SHA-256: `a4a28d28149156b9b81a0f75f0757291773d352f7752d6e23a356373e42ec9d7`.

## Frame cadence

**Pass aggregate**: peak frame render is **13.292 ms** (frame 920), mean live-frame render is **10.301 ms**, and **0/1503** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`sh_rasterize` averages 10.019 ms/f across complete post-setup windows; its worst window is 12.970 ms/f (frames 865–880).

Setup frame 1 rendered 17.967 ms. Published live-frame figures use complete counter windows through frame 1504; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:42) recorded peak 🟢 12.92 (24) and spilled 🟢 0/3487 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 913–928)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.437 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.143 | 1.00 | 0.23% |
| ↳ `sh_rasterize` | 12.968 | 1.00 | 20.77% |
| ↳ `sh_timeline_step` | 0.029 | 1.00 | 0.05% |
| ↳ `canvas_clear` | 0.083 | 1.00 | 0.13% |
| ↳ `canvas_buffer_wait` | 49.208 | 1.00 | 78.81% |

Wall min/avg/max = 62.202/62.456/62.558 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `sh_rasterize` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 21 | — | 13.292 | 0/64 | 12.968 | 3/4 |
| 20 | — | 13.281 | 0/64 | 12.970 | 3/4 |
| 22 | — | 13.269 | 0/64 | 11.843 | 3/4 |
| 19 | — | 12.136 | 0/64 | 11.820 | 3/4 |
| 23 | — | 12.064 | 0/64 | 10.501 | 3/4 |
| 13 | — | 11.623 | 0/64 | 11.260 | 3/4 |
| 12 | — | 11.619 | 0/64 | 11.296 | 3/4 |
| 14 | — | 11.532 | 0/64 | 10.052 | 3/4 |
| 18 | — | 10.822 | 0/64 | 10.506 | 3/4 |
| 24 | — | 10.760 | 0/64 | 9.811 | 3/4 |
| 11 | — | 10.365 | 0/64 | 10.014 | 3/4 |
| 15 | — | 10.253 | 0/64 | 9.339 | 3/4 |
| 17 | — | 10.127 | 0/64 | 9.803 | 3/4 |
| 1 | — | 10.067 | 0/64 | 9.013 | 3/4 |
| 16 | — | 9.826 | 0/64 | 9.487 | 3/4 |
| 7 | — | 9.823 | 0/62 | 9.476 | 2/3 |
| 6 | — | 9.817 | 0/33 | 9.449 | 2/2 |
| 8 | — | 9.781 | 0/64 | 8.883 | 3/4 |
| 10 | — | 9.671 | 0/64 | 9.351 | 3/4 |
| 9 | — | 9.365 | 0/64 | 9.049 | 3/4 |
| 2 | — | 9.209 | 0/64 | 8.325 | 3/4 |
| 5 | — | 9.205 | 0/64 | 8.868 | 3/4 |
| 4 | — | 8.895 | 0/64 | 8.580 | 3/4 |
| 3 | — | 8.676 | 0/64 | 8.340 | 3/4 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=999991 us
isr_wake       n=18433 cyc min/avg/max=302/1167/8158 ns min/avg/max=503/1945/13596 total=35876 us cpu=3.58%
isr_pack       n=2304 cyc min/avg/max=5244/5549/7045 ns min/avg/max=8740/9248/11741 total=21310 us cpu=2.13%
isr_dma_submit n=2304 cyc min/avg/max=464/560/607 ns min/avg/max=773/933/1011 total=2152 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `sh_rasterize` — 12.968 ms/f, 20.77% of the peak window.
2. `pov_preserve_half` — 0.143 ms/f, 0.23% of the peak window.
3. `canvas_clear` — 0.083 ms/f, 0.13% of the peak window.
4. `sh_timeline_step` — 0.029 ms/f, 0.05% of the peak window.

README cells: peak 🟢 13.292 (24), spilled 🟢 0/1503 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=SphericalHarmonics  windows=94
  [PASS] captured >=3 windows (94)
  [PASS] every window names one effect and resolution (SphericalHarmonics 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] cycle returns to its first index (peak idx 24, 24 distinct)
  [PASS] all 24 presets visited (24)
  [PASS] per-frame render is render, not wall (0 of 94 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 94 windows are short)
  [PASS] scope 'sh_rasterize' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (307.9 ppm, frames 865-880)

INVALID: SphericalHarmonics
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (307.9 ppm, frames 865-880). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=SphericalHarmonics`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
