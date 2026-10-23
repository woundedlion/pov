# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/hyperlattice_ship.log`, captured 2026-10-22 20:01 America/Los_Angeles on COM4.
Replaces `profile_hyperlattice_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 300 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh HyperLattice profile 300 16` |

Image size (`profile`, this effect only): `FLASH: code:143160, data:152396, headers:8572   free for files:1727488` / `RAM1: variables:312992, code:12888, padding:19880   free for local variables:178528` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/hyperlattice_ship_aa61e787859b_d73294c7627d`. Profile ELF SHA-256: `d73294c7627dff78488aa4c692f3350bde3c46068093a3513fac25d50148f39a`.

## Frame cadence

**Pass aggregate**: peak frame render is **55.080 ms** (frame 3209), mean live-frame render is **30.562 ms**, and **0/4767** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`hl_shader_draw` averages 28.113 ms/f across complete post-setup windows; its worst window is 44.094 ms/f (frames 3089–3104).

Setup frame 1 rendered 43.095 ms. Published live-frame figures use complete counter windows through frame 4768; a trailing partial window is excluded.

The previous shipping report (2026-10-06 20:43) recorded peak 🟢 51.77 (9) and spilled 🟢 0/5487 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 3201–3216)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.014 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `hl_shader_draw` | 43.537 | 1.00 | 70.21% |
| ↳ `hl_timeline_step` | 0.034 | 1.00 | 0.05% |
| ↳ `canvas_clear` | 0.098 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 16.074 | 1.00 | 25.92% |

Wall min/avg/max = 51.435/62.033/74.394 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `hl_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 6 | — | 55.080 | 0/559 | 44.094 | 34/35 |
| 3 | — | 39.201 | 0/559 | 29.101 | 34/35 |
| 9 | — | 37.673 | 0/416 | 33.969 | 25/26 |
| 8 | — | 34.224 | 0/679 | 28.852 | 42/42 |
| 4 | — | 33.832 | 0/439 | 28.725 | 26/27 |
| 5 | — | 33.002 | 0/679 | 29.971 | 42/43 |
| 2 | — | 30.248 | 0/679 | 25.721 | 42/43 |
| 7 | — | 30.234 | 0/439 | 25.107 | 26/27 |
| 1 | — | 27.017 | 0/318 | 23.814 | 18/19 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=994219 us
isr_wake       n=18327 cyc min/avg/max=372/1075/10257 ns min/avg/max=620/1791/17095 total=32853 us cpu=3.30%
isr_pack       n=2291 cyc min/avg/max=3594/4648/6680 ns min/avg/max=5990/7746/11133 total=17748 us cpu=1.78%
isr_dma_submit n=2291 cyc min/avg/max=463/559/618 ns min/avg/max=771/931/1030 total=2137 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `hl_shader_draw` — 43.537 ms/f, 70.21% of the peak window.
2. `pov_preserve_half` — 0.141 ms/f, 0.23% of the peak window.
3. `canvas_clear` — 0.098 ms/f, 0.16% of the peak window.
4. `hl_timeline_step` — 0.034 ms/f, 0.05% of the peak window.

README cells: peak 🟢 55.080 (9), spilled 🟢 0/4767 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=HyperLattice  windows=298
  [PASS] captured >=3 windows (298)
  [PASS] every window names one effect and resolution (HyperLattice 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 298 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 298 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): hl_shader_draw
  [PASS] scope 'hl_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (328.9 ppm, frames 3089-3104)

INVALID: HyperLattice
```

Validation limitations: cycle returns to its first index (peak idx 9, 8 distinct); all 9 presets visited (8); root cyc/600 vs wall sum within 5 ppm (328.9 ppm, frames 3089-3104). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
