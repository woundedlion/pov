# Shipping selective-O3 profiles

Ranked on-device results for the shipping `profile` image, covering the 34
effects in `HS_PHANTASM_EFFECT_LIST`. Peak is worst-frame render time; spilled
counts frames above the 62.5 ms display window. Colour is binary on the spilled
column: 🟢 zero spilled frames, 🔴 any spill. Cyclers (§) carry their preset
count after the peak.

| Effect | Dominant scope | Peak ms | Spilled | Captured |
|---|---|--:|--:|---|
| [DisplacementField](profile_displacementfield_teensy_2026-08-18.md) | fused ring-stack raster | 🟢 59.84 | 🟢 0/1408 (0%) | 2026-08-18 21:55 |
| [ShapeShifter](profile_shapeshifter_teensy_2026-08-08.md)§ | adaptive planar-star raster | 🟢 58.22 (9) | 🟢 0/2448 (0%) | 2026-08-08 17:54 |
| [HopfFibration](profile_hopffibration_teensy_2026-07-30.md) | trail raster + trail gate | 🟢 57.74 | 🟢 0/1088 (0%) | 2026-07-30 23:47 |
| [MeshFeedback](profile_meshfeedback_teensy_2026-08-05.md)§ | feedback flush (composite) | 🟢 57.70 (12) | 🟢 0/6688 (0%) | 2026-08-05 13:12 |
| [IslamicStars](profile_islamicstars_teensy_2026-07-28.md)§ | per-face SDF + opchain build legs | 🟢 56.91 (24) | 🟢 0/3328 (0%) | 2026-07-28 17:34 |
| [RingSpin](profile_ringspin_teensy_2026-07-25.md) | fused ring-group raster | 🟢 56.47 | 🟢 0/1088 (0%) | 2026-07-26 11:44 |
| [GSReactionDiffusion](profile_gsreactiondiffusion_teensy_2026-08-09.md) | integer opaque SSAA raster + sim | 🟢 56.28 | 🟢 0/2048 (0%) | 2026-08-09 16:34 |
| [Raymarch](profile_raymarch_teensy_2026-07-25.md) | volume ray-march | 🟢 52.99 | 🟢 0/1088 (0%) | 2026-07-26 11:38 |
| [BZReactionDiffusion](profile_bzreactiondiffusion_teensy_2026-08-03.md) | coefficient-factored SSAA raster | 🟢 50.70 | 🟢 0/2048 (0%) | 2026-08-03 00:33 |
| [VectorFacets](profile_vectorfacets_teensy_2026-08-16.md) ● | folded gnomonic dodecahedral vector mirror | 🟢 47.20 | 🟢 0/1088 (0%) | 2026-08-16 08:29 |
| [CurlLattice](profile_curllattice_teensy_2026-08-18.md)§ ● | curl-noise surface lattice | 🟢 45.18 (2) | 🟢 0/1728 (0%) | 2026-08-18 17:46 |
| [DreamBalls](profile_dreamballs_teensy_2026-08-09.md)§ | wireframe raster | 🟢 44.65 (5) | 🟢 0/3648 (0%) | 2026-08-09 18:37 |
| [HankinSolids](profile_hankinsolids_teensy_2026-07-25.md)§ | per-face SDF | 🟢 43.0 (19) | 🟢 0/3328 (0%) | 2026-07-26 11:55 |
| [Comets](profile_comets_teensy_2026-07-25.md)§ | point raster | 🟢 41.56 (12) | 🟢 0/4128 (0%) | 2026-07-26 11:43 |
| [MindSplatter](profile_mindsplatter_teensy_2026-08-07.md)§ | direct AA trail raster + clip gate | 🟢 38.95 (8) | 🟢 0/1728 (0%) | 2026-08-07 23:03 |
| [GnomonicStars](profile_gnomonicstars_teensy_2026-07-25.md) | star raster | 🟢 38.15 | 🟢 0/1088 (0%) | 2026-07-26 11:29 |
| [EquatorGrid](profile_equatorgrid_teensy_2026-08-17.md)§ ● | equirectangular dodecahedral grid | 🟢 36.51 (3) | 🟢 0/4128 (0%) | 2026-08-17 19:42 |
| [FacetWave](profile_facetwave_teensy_2026-08-16.md) ● | folded gnomonic dodecahedral wave mirror | 🟢 36.01 | 🟢 0/1088 (0%) | 2026-08-16 08:34 |
| [FacetGrid](profile_facetgrid_teensy_2026-08-16.md)§ ● | stereographic dodecahedral grid mirror | 🟢 35.52 (4) | 🟢 0/4128 (0%) | 2026-08-16 08:34 |
| [HexWave](profile_hexwave_teensy_2026-08-16.md) ● | stereographic hex-prism twin-wave | 🟢 34.45 | 🟢 0/1088 (0%) | 2026-08-16 08:36 |
| [SignalWeave](profile_signalweave_teensy_2026-08-16.md)§ ● | stereographic glitch wave-shear grid | 🟢 31.95 (4) | 🟢 0/4768 (0%) | 2026-08-16 08:44 |
| [PrismLattice](profile_prismlattice_teensy_2026-08-16.md) ● | stereographic prism polar lattice | 🟢 30.79 | 🟢 0/1088 (0%) | 2026-08-16 08:28 |
| [ContourLattice](profile_contourlattice_teensy_2026-08-16.md) ● | folded gnomonic affine lattice | 🟢 28.56 | 🟢 0/1088 (0%) | 2026-08-16 08:35 |
| [AlienOcean](profile_alienocean_teensy_2026-08-16.md) ● | folded gnomonic kaleidoscope grid | 🟢 28.20 | 🟢 0/1088 (0%) | 2026-08-16 08:30 |
| [KaleidoWave](profile_kaleidowave_teensy_2026-08-16.md) ● | stereographic kaleidoscope twin-wave | 🟢 27.89 | 🟢 0/1088 (0%) | 2026-08-16 08:28 |
| [CosmicEyeball](profile_cosmiceyeball_teensy_2026-08-16.md) ● | stereographic glitch mirror grid | 🟢 25.61 | 🟢 0/1088 (0%) | 2026-08-16 08:43 |
| [Fishbowl](profile_fishbowl_teensy_2026-08-02.md) | adaptive vertex build | 🟢 24.85 | 🟢 0/1088 (0%) | 2026-08-02 22:21 |
| [MobiusGrid](profile_mobiusgrid_teensy_2026-08-16.md)§ ● | stereographic Möbius twin-wave | 🟢 24.20 (2) | 🟢 0/2688 (0%) | 2026-08-16 08:46 |
| [GlitchGrid](profile_glitchgrid_teensy_2026-08-16.md) ● | folded gnomonic glitch mirror grid | 🟢 23.30 | 🟢 0/1088 (0%) | 2026-08-16 08:32 |
| [SphericalHarmonics](profile_sphericalharmonics_teensy_2026-08-21.md)§ | full-sphere harmonic shade | 🟢 12.29 (24) | 🟢 0/3488 (0%) | 2026-08-21 16:00 |
| [PetalFlow](profile_petalflow_teensy_2026-07-25.md) | ring raster | 🟢 11.71 | 🟢 0/1088 (0%) | 2026-07-26 11:37 |
| [Voronoi](profile_voronoi_teensy_2026-07-25.md) | block-union top-2 shade | 🟢 9.90 | 🟢 0/1088 (0%) | 2026-07-26 11:46 |
| [RingShower](profile_ringshower_teensy_2026-07-25.md) | ring raster | 🟢 4.07 | 🟢 0/1088 (0%) | 2026-07-26 11:40 |

**● captured 2026-08-16** — the fourteen composed effects first profiled in
that sweep.

**Eight rows carry a 2026-07-26 log** — RingSpin, Raymarch, HankinSolids,
Comets, GnomonicStars, PetalFlow, Voronoi and RingShower. Their per-effect
reports describe the 2026-07-25 sweep that preceded it, so each report's
headline peak is the earlier one and the peaks above are current.

§ ShapeShifter spans nine presets; its initial unlabeled frames and later
`Preset: 1/9` frames are one adaptive 208-count planar-star bucket.
MindSplatter spans eight presets. Each report folds its initial unlabeled
frames into preset 0. The five composed cyclers hold each preset for 600
frames and then morph over 480, so a preset owns 1,080 frames and their captures
are sized to wrap the full cycle.

This directory also holds the last captures of **Flyby** and **Liquid2D**, which
are not in the roster and are absent from the table above. ShaderWorkbench's captures
have been deleted: the composed-effect workbench migration (`69d4751c`) turned its
13-preset program bank into the fourteen ● effects.
