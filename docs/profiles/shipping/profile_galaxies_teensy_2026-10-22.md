# Galaxies on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/galaxies_ship.log`, captured 2026-10-22 18:55 America/Los_Angeles on COM4.
First canonical shipping report for Galaxies.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Galaxies 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 60 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh Galaxies profile 60 16` |

Image size (`profile`, this effect only): `FLASH: code:89272, data:152588, headers:9020   free for files:1780736` / `RAM1: variables:312960, code:23576, padding:9192   free for local variables:178560` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/galaxies_ship_aa61e787859b_e7c1159a7387`. Profile ELF SHA-256: `e7c1159a7387f719340f4060a196a627efdc7f7e691f347365503b033698029b`.

## Frame cadence

**Pass aggregate**: peak frame render is **44.932 ms** (frame 883), mean live-frame render is **24.663 ms**, and **0/927** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`gx_particle_step` averages 15.756 ms/f across complete post-setup windows; its worst window is 27.083 ms/f (frames 865–880).

Setup frame 1 rendered 1.521 ms. Published live-frame figures use complete counter windows through frame 928; a trailing partial window is excluded.

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 881–896)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.270 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.145 | 1.00 | 0.23% |
| ↳ `gx_draw_particles` | 15.941 | 1.00 | 25.60% |
| ↳ `gx_particle_step` | 27.071 | 1.00 | 43.47% |
| ↳ `gx_timeline_step` | 0.107 | 1.00 | 0.17% |
| ↳ `canvas_clear` | 0.096 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 18.863 | 1.00 | 30.29% |

Wall min/avg/max = 61.013/62.270/63.457 ms.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=998025 us
isr_wake       n=18397 cyc min/avg/max=348/1085/8770 ns min/avg/max=580/1808/14616 total=33284 us cpu=3.33%
isr_pack       n=2300 cyc min/avg/max=3687/4672/6193 ns min/avg/max=6145/7786/10321 total=17913 us cpu=1.79%
isr_dma_submit n=2300 cyc min/avg/max=411/565/2056 ns min/avg/max=685/941/3426 total=2167 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `gx_particle_step` — 27.071 ms/f, 43.47% of the peak window.
2. `gx_draw_particles` — 15.941 ms/f, 25.60% of the peak window.
3. `pov_preserve_half` — 0.145 ms/f, 0.23% of the peak window.
4. `gx_timeline_step` — 0.107 ms/f, 0.17% of the peak window.
5. `canvas_clear` — 0.096 ms/f, 0.15% of the peak window.

README cells: peak 🟢 44.932, spilled 🟢 0/927 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=Galaxies  windows=58
  [PASS] captured >=3 windows (58)
  [PASS] every window names one effect and resolution (Galaxies 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 58 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 58 windows are short)
  [PASS] scope 'gx_particle_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (11.4 ppm, frames 865-880)

INVALID: Galaxies
```

Validation limitations: root cyc/600 vs wall sum within 5 ppm (11.4 ppm, frames 865-880). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Galaxies`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
