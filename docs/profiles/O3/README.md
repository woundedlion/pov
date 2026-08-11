# Global-O3 reference profiles

Ranked paired captures for the single-effect `profile_o3` reference image.
Peak is worst-frame render time; spilled counts frames above the 62.5 ms
display window. Global O3 is a compiler ceiling, not a shippable roster image.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [ShaderBall](profile_shaderball_teensy_2026-08-11.md)§ ● | projection closure shader | 🔴 97.47 (5)<br>🟢 62.17 (24) | 🔴 168/371 (45.3%)<br>🟢 0/1837 (0%) | 2026-08-11 14:37 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-08-09.md) ● | integer opaque SSAA raster + sim | 🟢 56.97 | 🟢 0/2048 (0%) | 2026-08-09 16:37 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-08-08.md)§ ● | adaptive planar-star raster | 🟢 56.72 (9) | 🟢 0/2448 (0%) | 2026-08-08 17:57 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-08-03.md) ● | coefficient-factored SSAA raster | 🟢 50.90 | 🟢 0/2048 (0%) | 2026-08-03 00:36 |
| [DreamBalls](profile_dreamballs_teensy_2026-08-09.md)§ ● | wireframe raster | 🟢 42.94 (5) | 🟢 0/3648 (0%) | 2026-08-09 18:41 |
| [MindSplatter](profile_mindsplatter_teensy_2026-08-07.md)§ ● | direct AA trail raster + clip gate | 🟢 38.78 (8) | 🟢 0/1728 (0%) | 2026-08-07 23:02 |
| [ChaoticStrings](profile_chaoticstrings_teensy_2026-08-02.md) ● | adaptive vertex build | 🟢 22.16 | 🟢 0/1088 (0%) | 2026-08-02 22:23 |

§ ShapeShifter spans 10.08–56.72 ms across nine matched preset buckets; all
nine hold 16 fps.
The ShaderBall capture spans 29 presets and MindSplatter spans eight; each
report folds its initial unlabeled frames into preset 1.

**● refreshed 2026-08-11.**
