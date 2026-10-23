# KaleidoscopeFlowers on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopeflowers_ship.log`, captured 2026-10-22 20:19 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopeflowers_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeFlowers 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 173 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeFlowers profile 173 16` |

Image size (`profile`, this effect only): `FLASH: code:88568, data:155244, headers:9116   free for files:1778688` / `RAM1: variables:313024, code:17064, padding:15704   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/kaleidoscopeflowers_ship_aa61e787859b_2f65e9af94af`. Profile ELF SHA-256: `2f65e9af94afb231ecdbaff447e52db6c0997adb6f9cb78511b8343870499840`.

## Frame cadence

**Pass aggregate**: peak frame render is **33.982 ms** (frame 1957), mean live-frame render is **29.755 ms**, and **0/2735** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 26.220 ms/f across complete post-setup windows; its worst window is 29.282 ms/f (frames 1953–1968).

Setup frame 1 rendered 51.955 ms. Published live-frame figures use complete counter windows through frame 2736; a trailing partial window is excluded.

The previous shipping report (2026-10-06 19:04) recorded peak 🟢 31.42 (3) and spilled 🟢 0/4127 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 1953–1968)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.362 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.141 | 1.00 | 0.23% |
| ↳ `fx_shader_draw` | 29.282 | 1.00 | 46.96% |
| ↳ ↳ `fx_prepare_column` | 2.082 | 144.00 | 3.34% |
| ↳ `fx_prepare_frame` | 1.314 | 1.00 | 2.11% |
| ↳ `fx_advance` | 1.750 | 1.00 | 2.81% |
| ↳ `fx_timeline_step` | 0.195 | 1.00 | 0.31% |
| ↳ `canvas_clear` | 0.095 | 1.00 | 0.15% |
| ↳ `canvas_buffer_wait` | 29.227 | 1.00 | 46.87% |

Wall min/avg/max = 60.868/62.375/63.983 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 3 | — | 33.982 | 0/1058 | 29.282 | 66/66 |
| 2 | — | 32.317 | 0/1079 | 28.332 | 66/67 |
| 1 | — | 31.404 | 0/598 | 27.681 | 36/37 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=999811 us
isr_wake       n=18430 cyc min/avg/max=305/1207/16638 ns min/avg/max=508/2011/27730 total=37078 us cpu=3.70%
isr_pack       n=2304 cyc min/avg/max=5256/5677/7328 ns min/avg/max=8760/9461/12213 total=21803 us cpu=2.18%
isr_dma_submit n=2304 cyc min/avg/max=408/560/611 ns min/avg/max=680/933/1018 total=2152 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 29.282 ms/f, 46.96% of the peak window.
2. `fx_advance` — 1.750 ms/f, 2.81% of the peak window.
3. `fx_prepare_frame` — 1.314 ms/f, 2.11% of the peak window.
4. `fx_timeline_step` — 0.195 ms/f, 0.31% of the peak window.
5. `pov_preserve_half` — 0.141 ms/f, 0.23% of the peak window.
6. `canvas_clear` — 0.095 ms/f, 0.15% of the peak window.

README cells: peak 🟢 33.982 (3), spilled 🟢 0/2735 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=KaleidoscopeFlowers  windows=171
  [PASS] captured >=3 windows (171)
  [PASS] every window names one effect and resolution (KaleidoscopeFlowers 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 3, 2 distinct)
  [FAIL] all 3 presets visited (2)
  [PASS] per-frame render is render, not wall (0 of 171 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 171 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (213.3 ppm, frames 1953-1968)

INVALID: KaleidoscopeFlowers
```

Validation limitations: cycle returns to its first index (peak idx 3, 2 distinct); all 3 presets visited (2); root cyc/600 vs wall sum within 5 ppm (213.3 ppm, frames 1953-1968). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeFlowers`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
