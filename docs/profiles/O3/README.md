# Global-O3 reference profiles

Ranked paired captures for the single-effect `profile_o3` reference image.
Peak is worst-frame render time; spilled counts frames above the 62.5 ms
display window. Global O3 is a compiler reference, not a shippable roster
image or a guaranteed speed ceiling.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [DisplacementField](profile_displacementfield_teensy_2026-08-18.md) ● | fused ring-stack raster | 🟢 58.96 | 🟢 0/1408 (0%) | 2026-08-18 21:55 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-08-09.md) ● | integer opaque SSAA raster + sim | 🟢 56.97 | 🟢 0/2048 (0%) | 2026-08-09 16:37 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-08-08.md)§ ● | adaptive planar-star raster | 🟢 56.72 (9) | 🟢 0/2448 (0%) | 2026-08-08 17:57 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-08-03.md) ● | coefficient-factored SSAA raster | 🟢 50.90 | 🟢 0/2048 (0%) | 2026-08-03 00:36 |
| [AshCloud](profile_ashcloud_teensy_2026-08-23.md) ● | composed curl-noise shader | 🟢 48.35 | 🟢 0/1088 (0%) | 2026-08-23 22:06 |
| [HyperLattice](profile_hyperlattice_teensy_2026-08-25.md)§ | layered reflected-lattice shader | 🟢 47.79 (2) | 🟢 0/2688 (0%) | 2026-08-25 01:16 |
| [DreamBalls](profile_dreamballs_teensy_2026-08-09.md)§ ● | wireframe raster | 🟢 42.94 (5) | 🟢 0/3648 (0%) | 2026-08-09 18:41 |
| [MindSplatter](profile_mindsplatter_teensy_2026-08-07.md)§ ● | direct AA trail raster + clip gate | 🟢 38.78 (8) | 🟢 0/1728 (0%) | 2026-08-07 23:02 |
| [Fishbowl](profile_fishbowl_teensy_2026-08-02.md) ● | adaptive vertex build | 🟢 22.16 | 🟢 0/1088 (0%) | 2026-08-02 22:23 |

§ HyperLattice spans two presets. Both holds and both transition directions are
zero-spill; initial unlabeled frames fold into `cubic-flight`.
ShapeShifter spans 10.08–56.72 ms across nine
matched preset buckets; all nine hold 16 fps.
MindSplatter spans eight presets. Each report folds its initial unlabeled frames
into preset 0.

AshCloud is the first composed effect with a paired capture here. The other
fifteen composed effects have no O3-vs-shipping codegen pair on record.

**●** — capture refreshed 2026-08-15, except DisplacementField
(2026-08-18) and AshCloud (2026-08-23). Local to this index: the shipping
index marks composed effects with the same glyph.
