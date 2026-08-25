# On-device effect profiles — Teensy 4.0, segmented mode

On-device timing for the **36 effects in the Phantasm image**,
captured on bench-attached Teensy 4.0 boards running the shipping Phantasm
configuration (`POVSegmented<288, 4, 480>`, board = segment 0 master,
newlib-nano, DMA LEDs, flywheel + DMA ISRs live) via the `HS_PROFILE`
cycle-counter harness. `HS_EFFECT_LIST` holds 39 device effects, and
`HS_PHANTASM_EFFECT_LIST` excludes exactly three of them: Dynamo, MobiusRings
and Thrusters are Holosphere 96×20-only, never run in the 288×144 image, and
are not profiled here. The `Shader` authoring workbench and its `ShaderChain`
interpreter are simulator-only — `HS_ENABLE_SHADER_WORKBENCH` and
`HS_ENABLE_CHAIN_INTERPRETER` are never set for a device build — so they sit
outside that 39 rather than being excluded from it.

Each effect renders one **quadrant** ≈ **10,368 px**. A display window is
**62.5 ms**, so cadence quantizes: 16 fps (1 window), 8 fps (2), 5.3 fps (3).

Numbers are the shipping image only — the `-Os` `profile` env, whose `HS_O3`
regions activate on the `-Os` device build, so it measures the shipping
selective-O3 config by construction. `just profile <Effect>` regenerates one
row.

## How to read the table

**Peak** is the worst single frame's RENDER (frame minus the `*_buffer_wait`
display-sync idle). Wall time is render + sync wait, so its max is not a peak
render and is never used here. Per-frame averages are deliberately absent —
they hide the phase that sets cadence. The phase trees and per-preset spreads
live in each effect's report.

**Spilled** is `frames whose render overran one 62.5 ms window / total (%)`.

**Colour** is strict per phase: 🟢 zero spilled frames · 🟡 under 25% spilled ·
🔴 25% or more spilled. Peak ms is reported but does not set colour.

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
| [DisplacementField](shipping/profile_displacementfield_teensy_2026-08-18.md) | fused ring-stack raster | 🟢 59.84 | 🟢 0/1408 (0%) | 2026-08-18 21:55 |
| [ShapeShifter](shipping/profile_shapeshifter_teensy_2026-08-08.md)§ | adaptive planar-star raster | 🟢 58.22 (9) | 🟢 0/2448 (0%) | 2026-08-08 17:54 |
| [HopfFibration](shipping/profile_hopffibration_teensy_2026-07-30.md) | trail raster + trail gate | 🟢 57.74 | 🟢 0/1088 (0%) | 2026-07-30 23:47 |
| [MeshFeedback](shipping/profile_meshfeedback_teensy_2026-08-05.md)§ | feedback flush (composite) | 🟢 57.70 (12) | 🟢 0/6688 (0%) | 2026-08-05 13:12 |
| [IslamicStars](shipping/profile_islamicstars_teensy_2026-07-28.md)§ | per-face SDF + opchain build legs | 🟢 56.91 (24) | 🟢 0/3328 (0%) | 2026-07-28 17:34 |
| [RingSpin](shipping/profile_ringspin_teensy_2026-07-25.md) | fused ring-group raster (row-local walk) | 🟢 56.47 | 🟢 0/1088 (0%) | 2026-07-26 11:44 |
| [GSReactionDiffusion](shipping/profile_gsreactiondiffusion_teensy_2026-08-09.md) | integer opaque SSAA raster + sim | 🟢 56.28 | 🟢 0/2048 (0%) | 2026-08-09 16:34 |
| [Raymarch](shipping/profile_raymarch_teensy_2026-07-25.md) | volume ray-march (`-O3` march path) | 🟢 52.99 | 🟢 0/1088 (0%) | 2026-07-26 11:38 |
| [AshCloud](shipping/profile_ashcloud_teensy_2026-08-23.md) ● | composed curl-noise shader | 🟢 52.39 | 🟢 0/1088 (0%) | 2026-08-23 22:04 |
| [BZReactionDiffusion](shipping/profile_bzreactiondiffusion_teensy_2026-08-03.md) | coefficient-factored SSAA raster | 🟢 50.70 | 🟢 0/2048 (0%) | 2026-08-03 00:33 |
| [HyperLattice](shipping/profile_hyperlattice_teensy_2026-08-24.md)§ | layered reflected-lattice shader | 🟢 49.67 (2) | 🟢 0/2688 (0%) | 2026-08-25 01:12 |
| [KaleidoscopeStainedGlass](shipping/profile_kaleidoscopestainedglass_teensy_2026-08-16.md) ● | folded gnomonic dodecahedral vector mirror | 🟢 47.20 | 🟢 0/1088 (0%) | 2026-08-16 08:29 |
| [LatticeMelt](shipping/profile_latticemelt_teensy_2026-08-18.md)§ ● | curl-noise surface lattice | 🟢 45.18 (2) | 🟢 0/1728 (0%) | 2026-08-18 17:46 |
| [DreamBalls](shipping/profile_dreamballs_teensy_2026-08-09.md)§ | wireframe raster | 🟢 44.65 (5) | 🟢 0/3648 (0%) | 2026-08-09 18:37 |
| [HankinSolids](shipping/profile_hankinsolids_teensy_2026-07-25.md)§ | per-face SDF (`-O3` driver + `ConwayMorph`) | 🟢 43.0 (19) | 🟢 0/3328 (0%) | 2026-07-26 11:55 |
| [Comets](shipping/profile_comets_teensy_2026-07-25.md)§ | point raster | 🟢 41.56 (12) | 🟢 0/4128 (0%) | 2026-07-26 11:43 |
| [KaleidoscopeHexOil](shipping/profile_kaleidoscopehexoil_teensy_2026-08-17.md)§ ● | stereographic hex-prism spiral + direct noise | 🟢 40.52 (2) | 🟢 0/2208 (0%) | 2026-08-17 23:03 |
| [MindSplatter](shipping/profile_mindsplatter_teensy_2026-08-07.md)§ | direct AA trail raster + clip gate | 🟢 38.95 (8) | 🟢 0/1728 (0%) | 2026-08-07 23:03 |
| [GnomonicStars](shipping/profile_gnomonicstars_teensy_2026-07-25.md) | star raster | 🟢 38.15 | 🟢 0/1088 (0%) | 2026-07-26 11:29 |
| [KaleidoscopeFlowers](shipping/profile_kaleidoscopeflowers_teensy_2026-08-17.md)§ ● | equirectangular dodecahedral grid | 🟢 36.51 (3) | 🟢 0/4128 (0%) | 2026-08-17 19:42 |
| [KaleidoscopeMandala](shipping/profile_kaleidoscopemandala_teensy_2026-08-16.md) ● | folded gnomonic dodecahedral wave mirror | 🟢 36.01 | 🟢 0/1088 (0%) | 2026-08-16 08:34 |
| [KaleidoscopeSmooth](shipping/profile_kaleidoscopesmooth_teensy_2026-08-16.md)§ ● | stereographic dodecahedral grid mirror | 🟢 35.52 (4) | 🟢 0/4128 (0%) | 2026-08-16 08:34 |
| [KaleidoscopeHexBright](shipping/profile_kaleidoscopehexbright_teensy_2026-08-16.md) ● | stereographic hex-prism twin-wave | 🟢 34.45 | 🟢 0/1088 (0%) | 2026-08-16 08:36 |
| [AlienBrain](shipping/profile_alienbrain_teensy_2026-08-16.md)§ ● | stereographic glitch wave-shear grid | 🟢 31.95 (4) | 🟢 0/4768 (0%) | 2026-08-16 08:44 |
| [KaleidoscopePentBright](shipping/profile_kaleidoscopepentbright_teensy_2026-08-16.md) ● | stereographic prism polar lattice | 🟢 30.79 | 🟢 0/1088 (0%) | 2026-08-16 08:28 |
| [GridSpace](shipping/profile_gridspace_teensy_2026-08-16.md) ● | folded gnomonic affine lattice | 🟢 28.56 | 🟢 0/1088 (0%) | 2026-08-16 08:35 |
| [AlienOcean](shipping/profile_alienocean_teensy_2026-08-16.md) ● | folded gnomonic kaleidoscope grid | 🟢 28.20 | 🟢 0/1088 (0%) | 2026-08-16 08:30 |
| [KaleidoscopeHexSoft](shipping/profile_kaleidoscopehexsoft_teensy_2026-08-16.md) ● | stereographic kaleidoscope twin-wave | 🟢 27.89 | 🟢 0/1088 (0%) | 2026-08-16 08:28 |
| [CosmicEyeball](shipping/profile_cosmiceyeball_teensy_2026-08-16.md) ● | stereographic glitch mirror grid | 🟢 25.61 | 🟢 0/1088 (0%) | 2026-08-16 08:43 |
| [Fishbowl](shipping/profile_fishbowl_teensy_2026-08-02.md) | adaptive vertex build | 🟢 24.85 | 🟢 0/1088 (0%) | 2026-08-02 22:21 |
| [MobiusGrid](shipping/profile_mobiusgrid_teensy_2026-08-16.md)§ ● | stereographic Möbius twin-wave | 🟢 24.20 (2) | 🟢 0/2688 (0%) | 2026-08-16 08:46 |
| [AlienCore](shipping/profile_aliencore_teensy_2026-08-16.md) ● | folded gnomonic glitch mirror grid | 🟢 23.30 | 🟢 0/1088 (0%) | 2026-08-16 08:32 |
| [SphericalHarmonics](shipping/profile_sphericalharmonics_teensy_2026-08-21.md)§ | full-sphere harmonic shade | 🟢 12.29 (24) | 🟢 0/3488 (0%) | 2026-08-21 16:00 |
| [PetalFlow](shipping/profile_petalflow_teensy_2026-07-25.md) | ring raster | 🟢 11.71 | 🟢 0/1088 (0%) | 2026-07-26 11:37 |
| [Voronoi](shipping/profile_voronoi_teensy_2026-07-25.md) | block-union top-2 shade | 🟢 9.90 | 🟢 0/1088 (0%) | 2026-07-26 11:46 |
| [RingShower](shipping/profile_ringshower_teensy_2026-07-25.md) | ring raster | 🟢 4.07 | 🟢 0/1088 (0%) | 2026-07-26 11:40 |

**●** — the composed effects, in both tables on this page. Fourteen were
first profiled in the 2026-08-16 sweep; KaleidoscopeHexOil joined 2026-08-17 and
AshCloud 2026-08-23; LatticeMelt was re-captured 2026-08-18 as the
stage-families cut-over A/B gate. The other rows keep their own `Captured`
dates. The O3 index marks refresh dates with the same glyph; the two are
unrelated.

**Eight rows carry a 2026-07-26 log** — RingSpin, Raymarch, HankinSolids,
Comets, GnomonicStars, PetalFlow, Voronoi and RingShower. Their
per-effect reports describe the 2026-07-25 sweep that preceded it, so each
report's headline peak is the earlier one and the peaks above are current. `just
profile <Effect>` brings a report back level with its row.

## Paired shipping/O3 captures

Both peak columns precede the spill columns so the codegen delta reads
directly. Size deltas are O3 minus shipping.

Global O3 is a measurement reference, not a shipping candidate.

AshCloud is the first composed effect with a paired shipping/O3 capture. The
other fifteen composed effects have no codegen pair on record.

| Effect | Dominant scope | Ship peak ms | O3 peak ms | Ship spilled | O3 spilled | FLASH Δ | ITCM Δ | Captured |
|---|---|--:|--:|--:|--:|--:|--:|---|
| [DisplacementField](O3/profile_displacementfield_teensy_2026-08-18.md) | fused ring-stack raster | 🟢 59.84 | 🟢 58.96 | 🟢 0/1408 (0%) | 🟢 0/1408 (0%) | +25,216 B | +21,600 B | ship 2026-08-18 21:55<br>O3 2026-08-18 21:55 |
| [ShapeShifter](O3/profile_shapeshifter_teensy_2026-08-08.md)§ | adaptive planar-star raster | 🟢 58.22 (9) | 🟢 56.72 (9) | 🟢 0/2448 (0%) | 🟢 0/2448 (0%) | +28,616 B | +24,016 B | ship 2026-08-08 17:54<br>O3 2026-08-08 17:57 |
| [GSReactionDiffusion](O3/profile_gsreactiondiffusion_teensy_2026-08-09.md) | integer opaque SSAA raster + sim | 🟢 56.28 | 🟢 56.97 | 🟢 0/2048 (0%) | 🟢 0/2048 (0%) | +11,632 B | +10,624 B | ship 2026-08-09 16:34<br>O3 2026-08-09 16:37 |
| [AshCloud](O3/profile_ashcloud_teensy_2026-08-23.md) ● | composed curl-noise shader | 🟢 52.39 | 🟢 48.35 | 🟢 0/1088 (0%) | 🟢 0/1088 (0%) | +14,344 B | +11,888 B | ship 2026-08-23 22:04<br>O3 2026-08-23 22:06 |
| [BZReactionDiffusion](O3/profile_bzreactiondiffusion_teensy_2026-08-03.md) | coefficient-factored SSAA raster | 🟢 50.70 | 🟢 50.90 | 🟢 0/2048 (0%) | 🟢 0/2048 (0%) | +17,696 B | +16,256 B | ship 2026-08-03 00:33<br>O3 2026-08-03 00:36 |
| [HyperLattice](O3/profile_hyperlattice_teensy_2026-08-24.md)§ | layered reflected-lattice shader | 🟢 49.67 (2) | 🟢 47.79 (2) | 🟢 0/2688 (0%) | 🟢 0/2688 (0%) | +10,184 B | +7,936 B | ship 2026-08-25 01:12<br>O3 2026-08-25 01:16 |
| [DreamBalls](O3/profile_dreamballs_teensy_2026-08-09.md)§ | wireframe raster | 🟢 44.65 (5) | 🟢 42.94 (5) | 🟢 0/3648 (0%) | 🟢 0/3648 (0%) | +25,976 B | +16,272 B | ship 2026-08-09 18:37<br>O3 2026-08-09 18:41 |
| [MindSplatter](O3/profile_mindsplatter_teensy_2026-08-07.md)§ | direct AA trail raster + clip gate | 🟢 38.95 (8) | 🟢 38.78 (8) | 🟢 0/1728 (0%) | 🟢 0/1728 (0%) | +21,464 B | +18,832 B | ship 2026-08-07 23:03<br>O3 2026-08-07 23:02 |
| [Fishbowl](O3/profile_fishbowl_teensy_2026-08-02.md) | adaptive vertex build | 🟢 24.85 | 🟢 22.16 | 🟢 0/1088 (0%) | 🟢 0/1088 (0%) | +28,456 B | +20,688 B | ship 2026-08-02 22:21<br>O3 2026-08-02 22:23 |

**Nine pairs: six refreshed 2026-08-15, DisplacementField captured 2026-08-18,
AshCloud captured 2026-08-23, and HyperLattice completed 2026-08-25.**

## Captures of retired effects

[`retired/profile_flyby_teensy_2026-07-27.md`](retired/profile_flyby_teensy_2026-07-27.md)
and
[`retired/profile_liquid2d_teensy_2026-07-25.md`](retired/profile_liquid2d_teensy_2026-07-25.md)
preserve the last captures of Flyby and Liquid2D. Neither is in the roster, so
the current harness cannot regenerate them.

ShaderWorkbench's captures have been deleted. The composed-effect workbench migration
(`69d4751c`) turned its 13-preset program bank into the fourteen ● effects above
and left `Shader` as the authoring workbench rather than a shipping effect.

## Memory captures

[`memory/arena_high_water.md`](memory/arena_high_water.md) holds the per-effect
arena high-water marks against `DEVICE_GLOBAL_ARENA_SIZE`. Those come from a
host probe rather than the bench, so nothing in the timing tables above depends
on them.

## What the roster looks like

**Every current shipping capture is zero-spill.** HyperLattice's optimized 4D
transition now peaks at 49.67 ms with 0/2,688 spills, down from 65.37 ms and
104/2,592 spills. Global `-O3` peaks at 47.79 ms with zero spills and adds
10,184 B of flash code and 7,936 B of ITCM.

**The sixteen composed effects are green.** Their peaks run 23.30 ms
(AlienCore) to 52.39 ms (AshCloud), so the heaviest keeps 10.11 ms of the
62.5 ms display window.

The three slowest all perform per-sample noise work: AshCloud and LatticeMelt
integrate simplex curl-noise surfaces, while KaleidoscopeStainedGlass carries a
`Warp::VectorNoise` lookup in its outer warp.

Five of them cycle presets, and their spreads are narrow — KaleidoscopeFlowers 1.03×
across three presets, KaleidoscopeSmooth 1.19× across four, AlienBrain 1.03× across
four. A preset morph is a 480-frame parameter interpolation, not a structural
change, so a transition is bounded by its two endpoints rather than being a
cost peak of its own — a preset bank cannot introduce a cost outlier the way a
bank of separately compiled programs could.

AshCloud's global-O3 twin peaks at 48.35 ms, 7.7% below shipping while adding
14,344 B of flash code and 11,888 B of ITCM code. The other fifteen composed
effects have no paired global-O3 capture.

**MindSplatter completes all eight presets at 16 fps.** The shipping cycle has
0/1728 spills at a 38.95 ms peak; global O3 has the same zero-spill result at a
38.78 ms peak while adding 21,464 B of flash and 18,832 B of ITCM.

**MeshFeedback is green on all 12 styles** — 0/6688 at a 57.70 ms peak
(SlowDust), worst hold Smoke at 48.86 ms of flush. `feedback_composite` is 67%
of the frame. Span across styles is only 1.3×. `fbaf81e2` recovered 1.67 ms of a
regression `7846bc81` introduced; this metric sits on a 32-byte ITCM alignment
cliff, so treat sub-0.3 ms movements as code placement until a disassembly or a
same-image A/B says otherwise.

**GSReactionDiffusion now generates and rebakes a palette at every reaction
boundary.** Four captured rebakes cost 2.348–2.417 ms each. The shipping pass
peaks at 56.277 ms with 0/2048 spills, retaining 6.223 ms of display-window
margin, so cycling prebaked flash palettes is unnecessary.

**DisplacementField now sits just over budget** — a 59.84 ms peak against the
59.1 ms ISR-adjusted budget, about 1.3% above it, though it still spills
nothing (0/1408). Sixteen frames in two adjacent peak windows exceed the
62.5 ms display window. Global O3 brings the peak to 58.96 ms, back inside
the budget.

**IslamicStars is green on all 24 shapes** — 0/3328 at a 56.91 ms peak, with
the worst hold (`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`) at 47.05 ms and
the peak frame falling in a dual-bridge build leg. Its row was re-captured at
tip `542a5b49` (2026-07-28), 664 commits after the rest of the table's
`0bbc56e3`; the numbers moved under 0.1 ms.

Global O3 is 0.695 ms slower than GS shipping at 56.972 ms while adding
10,624 B of ITCM. The selective-O3 image is both smaller and faster, and its
generative-palette transition stays within the 16 fps budget.

**BZReactionDiffusion now peaks at 50.705 ms** with 0/2048 spills, down
8.378 ms (14.2%) from its matched pre-optimization capture. Its final
factorization is visually equivalent rather than framebuffer-bit-exact: a
256-frame production comparison measured 131.46 dB PSNR and no coverage
changes.

**ShapeShifter is green across all nine presets.** Its adaptive 208-contour
planar star reaches approximately 290-count density at the equator while
maintaining approximately 145-count density at the poles. The shipping image
peaks at 58.22 ms with 0/2448 spills; global O3 peaks at 56.72 ms and adds
28,616 B of flash code plus 24,016 B of ITCM code. Reusing the shaded fragment
shader removes 5,136 B from the full-roster ITCM image without a measured
shipping performance regression.
