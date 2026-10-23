# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-10-22, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/kaleidoscopehexoil_ship.log`, captured 2026-10-22 20:09 America/Los_Angeles on COM4.
Replaces `profile_kaleidoscopehexoil_teensy_2026-10-06.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, flywheel + DMA ISRs live, COM4 |
| Image | `profile`: shipping `-Os` base with selective `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, source `aa61e787859bf31d69057b436b08ad8a751f1bb4` |
| Method | `HS_PROFILE` cycle scopes, 16-frame windows, 105 s capture; current Phantasm show-segment budget |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-profiles HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 HS_PIO_JOBS=4 bash tools/profile_one.sh KaleidoscopeHexOil profile 105 16` |

Image size (`profile`, this effect only): `FLASH: code:88120, data:158796, headers:9084   free for files:1775616` / `RAM1: variables:313024, code:16856, padding:15912   free for local variables:178496` / `RAM2: variables:520064  free for malloc/new:4224`.

Capture provenance retains the source snapshot, compiler/ABI attestation, build logs and image artifacts in `build/prof/artifacts/kaleidoscopehexoil_ship_aa61e787859b_1c918ac86187`. Profile ELF SHA-256: `1c918ac861874d22abc1072d2901587f716533ed2efb95b1d8b18536cb64546f`.

## Frame cadence

**Pass aggregate**: peak frame render is **49.679 ms** (frame 891), mean live-frame render is **44.799 ms**, and **0/1647** frames spilled (0.00%). Setup frame 1 is excluded from those figures.
`fx_shader_draw` averages 41.167 ms/f across complete post-setup windows; its worst window is 44.765 ms/f (frames 769–784).

Setup frame 1 rendered 71.213 ms. Published live-frame figures use complete counter windows through frame 1648; a trailing partial window is excluded.

The previous shipping report (2026-10-06 18:47) recorded peak 🟢 37.08 (2) and spilled 🟢 0/2207 (0.00%).

A display window is 62.5 ms; render above that budget spills. `canvas_buffer_wait` is display-sync idle.

## Phase-by-phase readout

### Peak window (frames 881–896)

Scopes below are inclusive window averages. Nested rows overlap their parents; percentages use the root counter.

| Scope | ms/f | Calls/f | Root share |
|---|--:|--:|--:|
| `frame` | 62.364 | 1.00 | 100.00% |
| ↳ `pov_preserve_half` | 0.135 | 1.00 | 0.22% |
| ↳ `fx_shader_draw` | 44.512 | 1.00 | 71.37% |
| ↳ ↳ `fx_prepare_column` | 6.626 | 144.00 | 10.63% |
| ↳ `fx_prepare_frame` | 0.876 | 1.00 | 1.41% |
| ↳ `fx_advance` | 2.225 | 1.00 | 3.57% |
| ↳ `fx_timeline_step` | 0.209 | 1.00 | 0.34% |
| ↳ `canvas_clear` | 0.103 | 1.00 | 0.16% |
| ↳ `canvas_buffer_wait` | 14.153 | 1.00 | 22.69% |

Wall min/avg/max = 59.790/62.377/65.098 ms.

### Per-preset table

Parser-owned frame buckets include the transition following each entry; setup frame 1 is excluded. Clean holds use the modal call count, exclude windows straddling an advance, and omit the startup window. Startup may appear as a separate ownership bucket for named-shape markers.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `fx_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 49.679 | 0/1049 | 44.765 | 65/65 |
| 1 | — | 46.738 | 0/598 | 42.218 | 36/37 |

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot:

```text
ISR window=999936 us
isr_wake       n=18432 cyc min/avg/max=297/1183/15626 ns min/avg/max=495/1971/26043 total=36342 us cpu=3.63%
isr_pack       n=2304 cyc min/avg/max=4774/5511/7161 ns min/avg/max=7956/9185/11935 total=21165 us cpu=2.11%
isr_dma_submit n=2304 cyc min/avg/max=417/565/3939 ns min/avg/max=695/941/6565 total=2170 us cpu=0.21%
```

`isr_wake` includes packing and submission; their shares must not be added. CPU shares use the separately logged ISR interval. All cycle scopes include interrupt time because CYCCNT free-runs.

## Summary ranking

1. `fx_shader_draw` — 44.512 ms/f, 71.37% of the peak window.
2. `fx_advance` — 2.225 ms/f, 3.57% of the peak window.
3. `fx_prepare_frame` — 0.876 ms/f, 1.41% of the peak window.
4. `fx_timeline_step` — 0.209 ms/f, 0.34% of the peak window.
5. `pov_preserve_half` — 0.135 ms/f, 0.22% of the peak window.
6. `canvas_clear` — 0.103 ms/f, 0.16% of the peak window.

README cells: peak 🟢 49.679 (2), spilled 🟢 0/1647 (0.00%).

## Caveats

- The capture is limited to the effect’s current show segment; peaks describe the observed segment.
- Counter averages are inclusive and retain instrumentation overhead. MIXED-PARENT and DUPLICATE-NAME rows are not exclusive costs of the displayed parent.

Parser validation:

```text
effect=KaleidoscopeHexOil  windows=103
  [PASS] captured >=3 windows (103)
  [PASS] every window names one effect and resolution (KaleidoscopeHexOil 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 103 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 103 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (239.4 ppm, frames 769-784)

INVALID: KaleidoscopeHexOil
```

Validation limitations: cycle returns to its first index (peak idx 2, 1 distinct); all 2 presets visited (1); root cyc/600 vs wall sum within 5 ppm (239.4 ppm, frames 769-784). The raw capture is retained with these diagnostics; this report does not claim a fully validated cycle.

The ownership table includes the parser’s inferred startup entry. The validation visit count uses explicit advance markers and requires a subsequent return to the first index.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexOil`, `HS_PROFILE_WINDOW=16`; the Reproduce command builds, flashes and captures under the COM4 device lock.
