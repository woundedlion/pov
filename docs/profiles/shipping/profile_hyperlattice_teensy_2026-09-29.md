# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Raw capture: `build/prof/hyperlattice_ship.log` (standard cycle) and
`build/prof/hyperlattice_experimental_ship.log` (opt-in experimental cycle).
Replaces `profile_hyperlattice_teensy_2026-09-28.md`. Every per-change capture
of the optimization campaign is kept under
`build/prof/hyperlattice_campaign_2026-09-29/`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live; standard cycle on COM4, experimental cycle and pinned presets on COM3 (paired captures of one image on both boards agree within 0.06%) |
| Image | `profile` env; the cached shader uses `HS_HOT_FLASH_MEMBER` (`HS_O3_FN` plus `hot`), with no ITCM placement |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist; cycles captured at `b1ccbfa39`, pinned presets at the campaign branch before its rebase onto `67f935538` |
| Method | `HS_PROFILE` cycle scopes, window 16. Standard cycle: 120 s, `-D HS_PROFILE_EPOCH_REVS=1200`. Experimental cycle: 345 s, `-D HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -D HS_PROFILE_EPOCH_REVS=2800`. Per-preset A/B: 30 s pinned captures, `-D HS_PROFILE_PRESET=<i>` |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile 120 16 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size (standard profile image): `FLASH: code:77880, data:150016,
headers:8648` / `RAM1: variables:314944, code:13384, padding:19384, free:176576`
/ `RAM2: variables:520064, free:4224`.

Default Phantasm image against master `67f935538`: RAM1 code 168,856 →
168,936 B (+80 B, all inlining drift in other effects; HyperLattice's own ITCM
fell 240 B), FLASH code 512,264 → 529,344 B (+17,080 B). Experimental image:
RAM1 code 172,472 → 179,752 B (+7,280 B, again other effects' functions whose
inlining moved as the unit grew), FLASH code +88,720 B. RAM1 variables and the
12,896 B stack floor are unchanged in both.

Exactness cross-check: window frames 1185–1200 root counter cyc ÷ 600 MHz
matches the measured wall sum within **2.9 ppm**.

## Frame cadence

**Pass aggregate** (standard cycle): `hl_shader_draw` avg 23.76 ms/f, worst
window 28.09 ms/f (frames 1185–1200), peak frame render **36.65 ms** (frame
1198), spilled **0/1887**. Setup frame 1 is excluded.

A display window is 62.5 ms; the effect renders one quadrant, 146×73 = 10,658
rays per frame. Every frame of the standard cycle holds 16 fps with at least
25 ms of margin. The `canvas_buffer_wait` scope is the round-up idle to the
next display flip, by design.

## Phase-by-phase readout

Phase schedule: three presets, each a 320-frame hold followed by a 240-frame
segue. Cubic Flight and Cubic Wide Flight share a pattern and view, so that
segue interpolates every continuous control. The segues into and out of
Hypercube Flight change the view: each holds the outgoing preset while it dims
to black, switches in the dark, and brightens the incoming one.

### Fade into the hypercube (frames 1185–1200, worst window)

```
frame                  62.84 ms  37.70 Mcyc  100%
  hl_shader_draw       28.09 ms  16.85 Mcyc   45%  x1.0
  pov_preserve_half     0.14 ms   0.09 Mcyc    0%
  canvas_clear          0.09 ms   0.05 Mcyc    0%
  hl_timeline_step     15.6 us    9.4 kcyc     0%
  canvas_buffer_wait   32.76 ms  19.66 Mcyc   52%
```

Wall min/avg/max = 57.57/62.84/68.31 ms. The window sits in the second half of
the dip into Hypercube Flight, which renders the hypercube's own parameters;
its peak is a hold-cost camera position. Holds run 25–27 ms mean.

### Per-preset table

Pinned 30 s captures at the final code, frame 1 excluded, and the full
experimental cycle's buckets (captured 2026-09-29 16:11 on COM3, after presets
began carrying their own departures). A preset's marker fires when its
parameters are adopted: at the start of a morph, and at the dark midpoint of
a fade. So a bucket holds the morph into its preset, or the second half of the
fade into it, then its hold and the first half of any fade out of it.

| # | Preset | Pinned peak ms | Pinned mean ms | Cycle bucket peak ms | Cycle spilled |
|---|---|--:|--:|--:|--:|
| 5 | Octet 4D Flight | 56.12 | 47.00 | 🟢 55.87 | 0/559 |
| 8 | Shell 4D Flight | 51.17 | 35.04 | 🟢 36.76 | 0/559 |
| 2 | Hypercube Flight | 34.36 | 26.90 | 🟢 36.70 | 0/559 |
| 3 | Octet Flight | 31.92 | 28.12 | 🟢 32.23 | 0/439 |
| 7 | Shell Close Flight | 31.64 | 27.20 | 🟢 31.22 | 0/679 |
| 4 | Octet Wide Flight | 31.26 | 29.77 | 🟢 31.75 | 0/679 |
| 1 | Cubic Wide Flight | 29.15 | 25.59 | 🟢 29.64 | 0/817 |
| 6 | Shell Flight | 27.88 | 23.55 | 🟢 28.03 | 0/439 |
| 0 | Cubic Flight | 26.37 | 24.72 | 🟢 26.41 | 0/758 |

The experimental cycle wrapped to preset 1 (ten `Preset:` markers) and held
16 fps on every one of its 5,487 live frames. Octet 4D Flight's bucket peak,
55.87 ms at frame 3209, is that preset dimming at its own parameters; the 14:51
capture credited the same frame to Shell Flight, when the marker fired at the
start of the fade. Cubic Flight's bucket peak is setup frame 1; its peak past
setup is 26.41 ms. The wrap's morph into Cubic Wide Flight joins that preset's
bucket, hence its 817 frames. Root cycles match the wall sum within 0.6 ppm
(frames 3089–3104).

### Per-pixel figures

The shader writes premultiplied pixels directly; there is no `filter_blend`
counter. Standard cycle: 23.76 ms/f mean over 10,658 rays is about 1,338
cycles per ray.

## Column-ISR / DMA marshaling cost

```
isr_wake         1158/frame  min/avg/max 0.6/1.6/11.0 us  cpu 2.99%
isr_pack          145/frame  min/avg/max 6.2/6.7/9.2 us   cpu 1.54%
isr_dma_submit    145/frame  min/avg/max 0.7/0.9/1.0 us   cpu 0.21%
```

- Packing marshals the LED data on the CPU; submission only launches the DMA.
- A 600-byte transfer takes about 230 µs at the requested 24 MHz (LPSPI framing model), asynchronously.
- The ISRs take 4.7% of the CPU, leaving about 59.5 ms of render budget per
  62.5 ms window; neither cycle needs further speedup.

## Summary ranking

1. `hl_shader_draw` — 45% of the worst window, 28.09 ms/f.
2. `pov_preserve_half` — 0.2%, 0.14 ms/f.
3. `canvas_clear` — 0.1%, 0.09 ms/f.

Against the previous shipping report the standard cycle's peak fell from 56.12
to 36.65 ms, and both cycles hold 16 fps on every frame.

## Optimization ledger

Each change was profiled on device before landing. Peak ms are pinned 30 s
captures unless noted; "ITCM" is the profile image's RAM1 code delta, zero
throughout except where listed (the traced paths all run from cached flash).

| Change | Presets | Peak ms before → after | ITCM |
|---|---|---|--:|
| Shell layer march (3D, single owner) | Shell / Shell Close | 123.5 → 35.0 / 126.0 → 40.3 | 0 |
| Dedicated cubic merge loop | Hypercube / Cubic Wide / Cubic | 57.7 → 46.5 / 48.9 → 39.6 / 45.5 → 36.2 | 0 |
| Per-axis crossings, sort covered only, one division | same | 46.5 → 40.2 / 39.6 → 34.4 / 36.2 → 30.9 | 0 |
| Inline cubic scan lambda (reverted) | Hypercube / Cubic | ±0 | +12,880 |
| Octet per-stream walks, sort covered only | Octet / Octet Wide | 46.1 → 41.3 / 45.0 → 39.6 | 0 |
| Float palette composite | Hypercube / Octet / Shell Close | 40.2 → 37.0 / 41.3 → 39.2 / 40.3 → 36.0 | 0 |
| Direct shell scan | Shell / Shell Close | 35.0 → 29.4 / 36.0 → 33.3 | 0 |
| Direct 3D octet scan | Octet / Octet Wide | 39.2 → 36.6 / 39.6 → 35.9 | 0 |
| Fixed two-shell slice scan | Hypercube | 37.0 → 35.5 | 0 |
| Fixed two-shell cubic scan (with the palette composite) | Cubic Wide / Cubic | 34.4 → 30.0 / 30.9 → 27.1 | 0 |
| Static octet owners in registers | Octet / Octet Wide | 36.6 → 32.6 / 35.9 → 31.6 | 0 |
| Canonical-frame 4D octet | Octet 4D (old slice params) | 137.7 → 63.5 | 0 |
| Tighter 4D bound plus lazy parity (reverted) | Octet 4D (old params) | 63.5 → 67.0 | 0 |
| 4D shell neighbor march | Shell 4D | 179.8 → 59.9 | 0 |
| Stateless 4D octet plane bound | Octet 4D | 74.5 → 72.7 | 0 |
| Paired crossings (reverted) | Octet 4D | 72.7 → 76.6 | 0 |
| 16-bit fixed-point plane bound | Octet 4D | 72.7 → 67.8 | 0 |
| Packed halfword adds (reverted) | Octet 4D | 67.8 → 69.3 | 0 |
| Counted crossings, fixed-point class bound | Octet 4D | 67.8 → 65.7 | 0 |
| Across-residual bound | Octet 4D | 65.7 → 63.8 | 0 |
| Two-stage bound (reverted) | Octet 4D | 63.8 → 66.3 | 0 |
| Integer forward-differenced threshold (reverted) | Octet 4D | 63.8 → 67.5 | 0 |
| Shared family setup, one division | Octet 4D | 63.8 → 62.7 | 0 |
| Integer magnitude ranking | Octet 4D | 62.7 → 62.0 | 0 |
| Counted shell layers | Shell 4D / Shell Close / Shell | 59.8 → 56.6 / 33.5 → 32.2 / 29.5 → 28.5 | 0 |
| Lazy class bounds (reverted) | Octet 4D | 62.0 → 64.1 | 0 |
| `fminf` threshold clamp | Octet 4D | 62.0 → 61.0 | 0 |
| `fmaxf` shell offset | Shell 4D | 56.6 → 50.6 | 0 |
| `fminf`/`fmaxf` sweep of traced paths | all | e.g. Hypercube 35.5 → 34.4, Octet 4D 61.0 → 60.0; Shell 4D 50.6 → 51.2 | 0 |
| Branch-free class mask and across shift | Octet 4D | 60.0 → 56.1 | 0 |
| 32-bit fixed-point fractions (reverted) | Octet 4D | 56.1 → 61.4 | 0 |
| 3D neighbor march for wide lerped shells | Shell Close → Shell 4D segue | 142.7 → 51.0 (cycle bucket) | 0 |
| Midpoint switch for controls one side ignores (not landed) | Octet 4D → Shell segue | frozen pre-switch 68.0 → 64.6; cycle 65.5 → 76.7 (4D orientation moved) | 0 |
| Fade through black between pattern/view families | experimental cycle / standard cycle | 65.5 (9 spilled) → 56.0 (0 spilled) / 37.46 → 36.65 | 0 |

Octet 4D Flight and Shell 4D Flight joined the cycle with new parameters
during the campaign; each row measures the parameters current at its time.

## Caveats

- CYCCNT includes ISR time (4.7% of the CPU); it is part of every scope.
- No per-pixel scopes run in these captures; the campaign's deep captures
  (`HS_PROFILE_DEEP`) were attribution-only and are not timing figures.
- The traced shader runs from cached flash as `HS_HOT_FLASH_MEMBER` (`HS_O3_FN` plus `hot`), with no ITCM placement.
- The experimental presets require `HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1`;
  epoch stretches only lengthen the effect instance, never a frame's cost.
- The pinned captures ran the campaign branch before its rebase onto master,
  which added peer commits outside HyperLattice; both cycles ran `b1ccbfa39`.

## Harness

`targets/Profile/Profile.ino` knobs: `HS_PROFILE_TARGET=HyperLattice`,
`HS_PROFILE_WINDOW=16`, `HS_PROFILE_EPOCH_REVS`, `HS_PROFILE_PRESET`, and
`HS_ENABLE_HYPERLATTICE_EXPERIMENTS=1` for the experimental cycle. One-liner:
`just profile HyperLattice 120`.
