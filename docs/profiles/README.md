# On-device effect profiles — Teensy 4.0, segmented mode

On-device timing for the **21 effects in the Phantasm image**, captured on
bench-attached Teensy 4.0 boards running the shipping Phantasm configuration
(`POVSegmented<288, 4, 480>`, board = segment 0 master, newlib-nano, DMA LEDs,
flywheel + DMA ISRs live) via the `HS_PROFILE` cycle-counter harness. The full
effect roster is 23; Dynamo and Thrusters are Holosphere 96×20-only, never run
in the 288×144 Phantasm image (`HS_PHANTASM_EFFECT_LIST` excludes them), and are
not profiled here.

Each effect renders one **quadrant** ≈ **10,368 px**. A display window is
**62.5 ms**, so cadence quantizes: 16 fps (1 window), 8 fps (2), 5.3 fps (3).

**Full-roster sweep, 2026-07-26** (tip `0bbc56e3`), run in parallel across two
boards (COM3 + COM4). Numbers are the shipping image only — the `-Os` `profile`
env, whose `HS_O3` regions (docs/selective_o3_spec.md) activate on the `-Os`
device build, so it measures the shipping selective-O3 config by construction.
`just profile <Effect>` regenerates one row.

## How to read the table

**Peak** is the worst single frame's RENDER (frame minus the `*_buffer_wait`
display-sync idle). Wall time is render + sync wait, so its max is not a peak
render and is never used here. Per-frame averages are deliberately absent —
they hide the phase that sets cadence. The phase trees and per-preset spreads
live in each effect's report.

**Spilled** is `frames whose render overran one 62.5 ms window / total (%)`.

**Colour** is strict per phase: 🟢 zero spilled frames · 🟡 under 25% spilled ·
🔴 25% or more. Peak ms and spill fraction are both retained because cadence is
quantized.

**§ cyclers** carry one row per colour bucket, worst first, aligned across the
peak and spilled columns. A preset owns exactly the frames it was on screen for,
including its own transition in and out. These are **stricter** than the
clean-hold "N of M hold 16 fps" figure in the per-effect reports. The count in
parentheses after a bucket peak is presets in that bucket.

**Captured** is the raw log's mtime, local date and time. Rows are re-captured
piecemeal, so the roster never shares one timestamp — the date is how a reader
tells a fresh number from a stale one.

Rows rank by severity — worst colour bucket first, then spill fraction, then
peak.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [ShaderBall](shipping/profile_shaderball_teensy_2026-08-07.md)§ ● | stereographic closure shader | 🔴 97.93 (5)<br>🟢 61.23 (7) | 🔴 1501/2549 (59%)<br>🟢 0/2315 (0%) | 2026-08-07 23:12 |
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-07-28.md) | fused ring-stack raster | 🟢 58.71 | 🟢 0/1088 (0%) | 2026-07-28 17:41 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-07-30.md) | trail raster + trail gate | 🟢 57.74 | 🟢 0/1088 (0%) | 2026-07-30 23:47 |
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-08-05.md)§ | feedback flush (composite) | 🟢 57.70 (12) | 🟢 0/6688 (0%) | 2026-08-05 13:12 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-08-04.md)§ ● | adaptive planar-star raster | 🟢 57.66 (9) | 🟢 0/2368 (0%) | 2026-08-04 11:06 |
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-07-28.md)§ | per-face SDF + opchain build legs | 🟢 56.91 (24) | 🟢 0/3328 (0%) | 2026-07-28 17:34 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-07-25.md) | fused ring-group raster (row-local walk) | 🟢 56.47 | 🟢 0/1088 (0%) | 2026-07-26 11:44 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-08-03.md) ● | integer opaque SSAA raster + sim | 🟢 55.71 | 🟢 0/2048 (0%) | 2026-08-03 00:27 |
| [Raymarch](shipping/profile_raymarch_teensy_2026-07-25.md) | volume ray-march (`-O3` march path) | 🟢 52.99 | 🟢 0/1088 (0%) | 2026-07-26 11:38 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-08-03.md) ● | coefficient-factored SSAA raster | 🟢 50.70 | 🟢 0/2048 (0%) | 2026-08-03 00:33 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-08-05.md)§ | wireframe raster | 🟢 44.88 (5) | 🟢 0/3408 (0%) | 2026-08-05 13:10 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-07-25.md)§ | per-face SDF (`-O3` driver + `ConwayMorph`) | 🟢 43.0 (19) | 🟢 0/3328 (0%) | 2026-07-26 11:55 |
| [Comets](shipping/profile_comets_teensy_2026-07-25.md)§ | point raster | 🟢 41.56 (12) | 🟢 0/4128 (0%) | 2026-07-26 11:43 |
| [MindSplatter](shipping/profile_mindsplatter_teensy_2026-08-07.md)§ ● | direct AA trail raster + clip gate | 🟢 38.95 (8) | 🟢 0/1728 (0%) | 2026-08-07 23:03 |
| [GnomonicStars](shipping/profile_gnomonicstars_teensy_2026-07-25.md) | star raster | 🟢 38.15 | 🟢 0/1088 (0%) | 2026-07-26 11:29 |
| [ChaoticStrings](shipping/profile_chaoticstrings_teensy_2026-08-02.md) ● | adaptive vertex build | 🟢 24.85 | 🟢 0/1088 (0%) | 2026-08-02 22:21 |
| [MobiusGrid](shipping/profile_mobiusgrid_teensy_2026-07-25.md) | curve raster | 🟢 17.24 | 🟢 0/1088 (0%) | 2026-07-26 11:35 |
| [SphericalHarmonics](shipping/profile_sphericalharmonics_teensy_2026-07-25.md)§ | field raster | 🟢 15.9 (24) | 🟢 0/3488 (0%) | 2026-07-26 11:59 |
| [PetalFlow](shipping/profile_petalflow_teensy_2026-07-25.md) | ring raster | 🟢 11.71 | 🟢 0/1088 (0%) | 2026-07-26 11:37 |
| [Voronoi](shipping/profile_voronoi_teensy_2026-07-25.md) | block-union top-2 shade | 🟢 9.90 | 🟢 0/1088 (0%) | 2026-07-26 11:46 |
| [RingShower](shipping/profile_ringshower_teensy_2026-07-25.md) | ring raster | 🟢 4.07 | 🟢 0/1088 (0%) | 2026-07-26 11:40 |

## Paired shipping/O3 captures

Both peak columns precede the spill columns so the codegen delta reads
directly. Size deltas are O3 minus shipping.

| Effect | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured |
|---|---|--:|--:|--:|--:|--:|--:|---|
| [ShaderBall](O3/profile_shaderball_teensy_2026-08-07.md)§ ● | stereographic closure shader | 🔴 97.93 (5)<br>🟢 61.23 (7) | 🔴 89.21 (3)<br>🟢 56.69 (9) | 🔴 1501/2549 (59%)<br>🟢 0/2315 (0%) | 🔴 1437/1593 (90%)<br>🟢 0/3335 (0%) | +29,376 B | +26,592 B | ship 2026-08-07 23:12<br>O3 2026-08-07 23:22 |
| [ShapeShifter](O3/profile_shapeshifter_teensy_2026-08-04.md)§ ● | adaptive planar-star raster | 🟢 57.66 (9) | 🟢 59.47 (9) | 🟢 0/2368 (0%) | 🟢 0/2368 (0%) | +30,088 B | +25,456 B | ship 2026-08-04 11:06<br>O3 2026-08-04 11:09 |
| [GSReactionDiffusion](O3/profile_gsreactiondiffusion_teensy_2026-08-03.md) ● | integer opaque SSAA raster + sim | 🟢 55.71 | 🟢 56.68 | 🟢 0/2048 (0%) | 🟢 0/2048 (0%) | +17,528 B | +16,432 B | ship 2026-08-03 00:27<br>O3 2026-08-03 00:31 |
| [BZReactionDiffusion](O3/profile_bzreactiondiffusion_teensy_2026-08-03.md) ● | coefficient-factored SSAA raster | 🟢 50.70 | 🟢 50.90 | 🟢 0/2048 (0%) | 🟢 0/2048 (0%) | +17,696 B | +16,256 B | ship 2026-08-03 00:33<br>O3 2026-08-03 00:36 |
| [MindSplatter](O3/profile_mindsplatter_teensy_2026-08-07.md)§ ● | direct AA trail raster + clip gate | 🟢 38.95 (8) | 🟢 38.78 (8) | 🟢 0/1728 (0%) | 🟢 0/1728 (0%) | +21,464 B | +18,832 B | ship 2026-08-07 23:03<br>O3 2026-08-07 23:02 |
| [ChaoticStrings](O3/profile_chaoticstrings_teensy_2026-08-02.md) ● | adaptive vertex build | 🟢 24.85 | 🟢 22.16 | 🟢 0/1088 (0%) | 🟢 0/1088 (0%) | +28,456 B | +20,688 B | ship 2026-08-02 22:21<br>O3 2026-08-02 22:23 |

**● refreshed 2026-08-07.**

## What the roster looks like

**Twenty of the twenty-one effects spill nothing** in their current shipping
captures. ShaderBall is the exception: five of its 12 shipping presets spill,
with 1501/2549 frames (59%) in the red bucket.

**MindSplatter completes all eight presets at 16 fps.** The shipping cycle has
0/1728 spills at a 38.95 ms peak; global O3 has the same zero-spill result at a
38.78 ms peak while adding 21,464 B of flash and 18,832 B of ITCM.

**ShaderBall exposes the merged 12-preset stereographic workload.** Shipping
peaks at 97.93 ms with five red presets. Global O3 lowers the peak to 89.21 ms
and leaves three red presets, at a cost of 29,376 B of flash and 26,592 B of
ITCM.

**MeshFeedback is green on all 12 styles** — 0/6688 at a 57.70 ms peak
(SlowDust), worst hold Smoke at 48.86 ms of flush. `feedback_composite` is 67%
of the frame. Span across styles is only 1.3×. `fbaf81e2` recovered 1.67 ms of a
regression `7846bc81` introduced; this metric sits on a 32-byte ITCM alignment
cliff, so treat sub-0.3 ms movements as code placement until a disassembly or a
same-image A/B says otherwise.

**GSReactionDiffusion now peaks at 55.708 ms** with 0/2048 spills, down from
94.830 ms and 1011/1024 spills. Typed rendering, one shared stencil per pixel,
safe seed acceptance, selective-O3 physics, six scaled integrations, and exact
integer opaque accumulation cut the peak by 41.3%. It is 0.708 ms above the
55 ms stretch target while retaining 6.792 ms of display-window margin.

**DisplacementField remains close behind** — 58.71 ms peak against a 59.4 ms
ISR-adjusted budget. It improved from 60.59 ms since 2026-07-25
(`filter_blend` fell 169 → 51 cyc/blend), but a ~1% workload increase would
still put it over.

**IslamicStars is green on all 24 shapes** — 0/3328 at a 56.91 ms peak, with
the worst hold (`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`) at 47.05 ms and
the peak frame falling in a dual-bridge build leg. Its row was re-captured at
tip `542a5b49` (2026-07-28), 664 commits after the rest of the table's
`0bbc56e3`; the numbers moved under 0.1 ms.

Global O3 is 0.969 ms slower than GS shipping at 56.677 ms while adding
16,432 B of ITCM. The selective-O3 image is both smaller and faster for the
final opaque-accumulation path.

**BZReactionDiffusion now peaks at 50.705 ms** with 0/2048 spills, down
8.378 ms (14.2%) from its matched pre-optimization capture. Its final
factorization is visually equivalent rather than framebuffer-bit-exact: a
256-frame production comparison measured 131.46 dB PSNR and no coverage
changes.

**ShapeShifter is green across all nine presets.** Its adaptive 208-contour
planar star reaches approximately 290-count density at the equator while
maintaining approximately 145-count density at the poles. The shipping image
peaks at 57.66 ms with 0/2368 spills; global O3 peaks at 59.47 ms and adds
30,088 B of flash code plus 25,456 B of ITCM code.
