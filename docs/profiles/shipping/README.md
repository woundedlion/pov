# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image. Peak is worst-frame
render time; spilled counts frames above the 62.5 ms display window. Cyclers
(§) show the worst cadence bucket and its preset count.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [ShaderBall](profile_shaderball_teensy_2026-08-07.md)§ ● | stereographic closure shader | 🔴 97.93 (5)<br>🟢 61.23 (7) | 🔴 1501/2549 (59%)<br>🟢 0/2315 (0%) | 2026-08-07 23:12 |
| [DisplacementField](profile_displacementfield_teensy_2026-07-28.md) | fused ring-stack raster | 🟢 58.71 | 🟢 0/1088 (0%) | 2026-07-28 17:41 |
| [HopfFibration](profile_hopffibration_teensy_2026-07-30.md) | trail raster + trail gate | 🟢 57.74 | 🟢 0/1088 (0%) | 2026-07-30 23:47 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-08-04.md)§ ● | adaptive planar-star raster | 🟢 57.66 (9) | 🟢 0/2368 (0%) | 2026-08-04 11:06 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-08-05.md)§ | feedback flush (composite) | 🟢 57.70 (12) | 🟢 0/6688 (0%) | 2026-08-05 13:12 |
| [IslamicStars](profile_islamicstars_teensy_2026-07-28.md)§ | per-face SDF + opchain build legs | 🟢 56.91 (24) | 🟢 0/3328 (0%) | 2026-07-28 17:34 |
| [RingSpin](profile_ringspin_teensy_2026-07-25.md) | fused ring-group raster | 🟢 56.47 | 🟢 0/1088 (0%) | 2026-07-26 11:44 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-08-03.md) ● | integer opaque SSAA raster + sim | 🟢 55.71 | 🟢 0/2048 (0%) | 2026-08-03 00:27 |
| [DreamBalls](profile_dreamballs_teensy_2026-08-05.md)§ | wireframe raster | 🟢 44.88 (5) | 🟢 0/3408 (0%) | 2026-08-05 13:10 |
| [Raymarch](profile_raymarch_teensy_2026-07-25.md) | volume ray-march | 🟢 52.99 | 🟢 0/1088 (0%) | 2026-07-26 11:38 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-08-03.md) ● | coefficient-factored SSAA raster | 🟢 50.70 | 🟢 0/2048 (0%) | 2026-08-03 00:33 |
| [HankinSolids](profile_hankinsolids_teensy_2026-07-25.md)§ | per-face SDF | 🟢 43.0 (19) | 🟢 0/3328 (0%) | 2026-07-26 11:55 |
| [Comets](profile_comets_teensy_2026-07-25.md)§ | point raster | 🟢 41.56 (12) | 🟢 0/4128 (0%) | 2026-07-26 11:43 |
| [MindSplatter](profile_mindsplatter_teensy_2026-08-07.md)§ ● | direct AA trail raster + clip gate | 🟢 38.95 (8) | 🟢 0/1728 (0%) | 2026-08-07 23:03 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-07-25.md) | star raster | 🟢 38.15 | 🟢 0/1088 (0%) | 2026-07-26 11:29 |
| [ChaoticStrings](profile_chaoticstrings_teensy_2026-08-02.md) ● | adaptive vertex build | 🟢 24.85 | 🟢 0/1088 (0%) | 2026-08-02 22:21 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-07-25.md) | curve raster | 🟢 17.24 | 🟢 0/1088 (0%) | 2026-07-26 11:35 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-07-25.md)§ | field raster | 🟢 15.9 (24) | 🟢 0/3488 (0%) | 2026-07-26 11:59 |
| [PetalFlow](profile_petalflow_teensy_2026-07-25.md) | ring raster | 🟢 11.71 | 🟢 0/1088 (0%) | 2026-07-26 11:37 |
| [Voronoi](profile_voronoi_teensy_2026-07-25.md) | block-union top-2 shade | 🟢 9.90 | 🟢 0/1088 (0%) | 2026-07-26 11:46 |
| [RingShower](profile_ringshower_teensy_2026-07-25.md) | ring raster | 🟢 4.07 | 🟢 0/1088 (0%) | 2026-07-26 11:40 |

§ ShapeShifter spans nine presets; its initial unlabeled frames and later
`Preset: 1/9` frames are one adaptive 208-count planar-star bucket.
ShaderBall spans 12 presets and MindSplatter spans eight; each report folds its
initial unlabeled frames into preset 1.

**● refreshed 2026-08-07.**
