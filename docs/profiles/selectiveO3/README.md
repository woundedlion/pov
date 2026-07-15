# Selective -O3 on-device profiles — Teensy 4.0, segmented mode (2026-07-15)

The shipping Phantasm image is `-Os`; the selective -O3 mechanism
(`HS_O3_BEGIN`/`HS_O3_END`/`HS_O3_FN`, `core/engine/platform.h`,
docs/selective_o3_spec.md) compiles only the measured hot loops at
`-O3 -ffast-math -fno-finite-math-only` while the rest of the image stays
`-Os`. Landed region set (commits `683d0eda`, `0f7b0616`, `e98f4651`):

- the terminal `Pipeline<W,H>` blend sink (`filter_blend` leaf),
- the fused scan drivers `Scan::RingGroup` (RingSpin) and
  `Scan::DistortedRingStack` (DisplacementField),
- `SDF::DistortedRing::distance_from_frame` / `polyline_distance`,
- `DisplacementField::draw_rings` (bake + hue prep + shader lambdas).

Cost on the shipping image: ITCM code 149,480 → **156,472 B** (+6,992 B),
intra-bank padding 7,368 B left (≥ the 4,096 B growth reserve), gated by the
`ram1.components.code` ceiling. Every other build (holosphere -O3, host,
WASM) is unaffected — the macros are no-ops there.

## How close does selective -O3 get to the -O3 ceiling?

All three configurations were captured on the same code state (pre-region
tip `2c2470b2` for `-Os` and the global-`-O3` twin; region tip for
selective), same device, same session. Hot-scope render time per frame:

| Effect / regime | Hot scope | -Os | **selective -O3** | -O3 ceiling | gap closed |
|---|---|--:|--:|--:|--:|
| RingSpin, median window | `rs_ring_scan` | 46.3 ms | **34.6 ms** | 33.3 ms | **90%** |
| RingSpin, heaviest window | `rs_ring_scan` | 57.5 ms | **41.0 ms** | 39.2 ms | **90%** |
| DisplacementField, NOISE-dwell median | `df_fused_scan` | 58.5 ms | **45.7 ms** | 44.0 ms | **88%** |
| shared blend leaf (RingSpin) | `filter_blend` | 149 cyc/blend | **84 cyc/blend** | 103 cyc/blend | >100% |

Cadence — the number that matters (62.5 ms display windows):

- **RingSpin: 8↔16 fps jitter → 16 fps hard-locked in every captured window**,
  identical to the global--O3 ceiling, with 33–43% of the frame back in
  `rs_buffer_wait` idle.
- **DisplacementField: 67% of the cycle now locks 16 fps** (BALLS phase,
  NOISE ramps, and the steady NOISE regime; all were 8 fps at `-Os` through
  the dwell). The **peak** NOISE-dwell stretch (~14 of 55 windows) still runs
  2 windows (8 fps): its render is 63–66 ms vs the window's 62.5. The global
  -O3 ceiling is itself marginal there — 59.7 ms render, window average
  62.56 ms with individual frames spilling to 65 ms — so the last stretch is
  not reachable by optimization level alone; it needs the render itself
  trimmed (~2–4 ms: the remaining -Os OKLab hue chain and
  `FastNoiseLite`-driven bake are the candidates, both measured; a hue-member
  `HS_O3_FN` promotion without the chain was measured dead and reverted,
  `7284bb09`).

Full standard-format reports:

- [profile_ringspin_teensy_2026-07-15.md](profile_ringspin_teensy_2026-07-15.md)
- [profile_displacementfield_teensy_2026-07-15.md](profile_displacementfield_teensy_2026-07-15.md)

Raw logs: `build/profile_ringspin_selo3.log`,
`build/profile_displacementfield_selo3_v2.log` (plus `_base_os` /
`_base_o3` twins for the same-tip baselines).

Point-in-time snapshots: `just profile <Effect>` now measures the shipping
selective-O3 configuration by construction (the `profile` env is `-Os`, so
the regions are active); `profile_o3` remains the untouched global-O3
reference.
