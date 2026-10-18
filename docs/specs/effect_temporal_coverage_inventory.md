# Effect temporal coverage inventory

Status: AUDITED on 2026-10-17 against source 7bf483103. Inventory only;
implementation choices and rollout order remain open.

## Result

The authored roster contains **43 effects**: **18 have column-time rendering
and 25 do not**. The optional ShaderChain interpreter adds a 44th path without
column-time rendering. The roster source is [targets/effects.h](../../targets/effects.h).

DreamBalls and every other animated forward renderer belong in the outstanding
inventory. They advance motion, then draw geometry from that frame's state;
they do not evaluate it at each destination column's emission time. This applies
to curves and particles as well as filled meshes. A primitive implemented using
a per-pixel distance function is still temporally uncovered if its transform
and shape state are held for the whole draw.

The same gap affects frame-wide inverse shaders, evolving simulation fields,
and feedback images. Spatial antialiasing, orientation history, interpolated
trail segments, and motion blur do not establish column-time continuity.

## Reading the inventory

**Covered** means the effect opts into the implemented composed temporal path:
geometry, continuous controls, clocks and color resources are sampled using
absolute column time. This is source-verified coverage, not a guarantee that
every spatial artifact or discrete transition is absent.

**Uncovered** means the rendering path lacks that timing mechanism. During
animation, samples on opposite sides of a sweep boundary can therefore have a
frame-sized state change instead of the ordinary column-sized change. The
existence of the gap is established from source; its visible severity depends
on motion, contrast, preset and observer. DreamBalls was specifically identified
by the user. The other entries have not been individually ranked by a new
visual/device survey. They are outstanding effects, not presumed seam-free
until someone reports them.

The state column records what currently changes or is resolved per frame. It
does not prescribe an interpolation policy. The linked effect file contains
the draw entry point and its animation/render callbacks.

## Uncovered: forward meshes

| Effect | Current rendering and changing state | Existing history or related machinery |
|---|---|---|
| [DreamBalls](../../effects/DreamBalls.h) | Orbiting wireframe copies through Plot::Mesh. Global orientation, orbit phase, tangent displacement, weave framing/gaps, per-sprite controls and crossfade opacity are resolved for a frame. Both original and woven/medial paths are uncovered. | Two sprite parameter/palette slots protect overlapping sprites; they are not consecutive temporal render endpoints. Rendering occurs inside timeline callbacks. |
| [HankinSolids](../../effects/HankinSolids.h) | Filled face/strap scans through Scan::Mesh. Camera orientation, Hankin shape evolution, Conway morph geometry, strap opening/closing and per-face palette handoffs change. | Compiled morph legs and transition palettes preserve authored shape changes; one rotated mesh and shading state supply each draw. |
| [IslamicStars](../../effects/IslamicStars.h) | Filled polyhedral scans through Scan::Mesh. Rotation, ripple displacement, recipe-build/morph geometry, opacity and face shading change. | Ripple preparation and compiled build legs are frame-based; timeline callbacks render held shapes and build frames. |

## Uncovered: forward rings, contours and stars

| Effect | Current rendering and changing state | Existing history or related machinery |
|---|---|---|
| [DisplacementField](../../effects/DisplacementField.h) | Distorted ring stack. Noise phase, falling-ball displacement, color spin, field gains and ball lifetimes change before frame-wide ring bakes/scanning. | Cached knots, chunks and prepared displacement fields accelerate the current frame. |
| [GnomonicStars](../../effects/GnomonicStars.h) | Scattered polygon-star scans. Animated warp and orientation move the stars; live radius, side and point controls alter their geometry. | Base Fibonacci positions are cached; transformed positions and bases are held for each star draw. |
| [MobiusRings](../../effects/MobiusRings.h) | Rasterized latitude/longitude curves. Mobius warp, global rotation, counter-rotation, moving pole holes, scrolling color coordinates and palette wipe change. | Sampled curve fragments and pole guards address spatial geometry, not emission time. |
| [PetalFlow](../../effects/PetalFlow.h) | Rasterized petal-shaped rings. Log-radial ring positions, orientation, spawn/retirement, opacity and palette coordinates evolve. | Stable live-ring records and cached petal wobble are used by an update-and-draw frame pass. |
| [RingShower](../../effects/RingShower.h) | Rasterized expanding rings. Each slot's age controls radius and fade-in; ring spawning, expiration and palette-slot reuse change the scene. | Per-ring normal, lifetime and baked palette persist, but age is sampled once for its draw. |
| [ShapeShifter](../../effects/ShapeShifter.h) | Concentric circle/star/shape scans. Waveform phase, radial sweep, contour radii, distribution, opacity and live controls change. | Adaptive contour/scan paths use the same advanced frame phase. |
| [Thrusters](../../effects/Thrusters.h) | Distorted main ring and spawned opposing rings. Global ring waveform, rotation/warp, fire events and each thruster's growth/fade age change. | Spawned contexts retain orientation and integer age; the modulo waveform clock and ages advance once per frame. |

## Uncovered: trails and replicated strands

| Effect | Current rendering and changing state | Existing history or related machinery |
|---|---|---|
| [Comets](../../effects/Comets.h) | Scanned points sampled along an orientation trail. The Lissajous head, global orientation, trail placement/fade and palette wipe change. | Orientation history and deep tweening construct the trail, then all its samples are drawn from frame-wide state. |
| [Dynamo](../../effects/Dynamo.h) | Replicated trailing strand nodes through world trail/replication/orientation filters. Strand pulls, emission count, orientation, trail ages and palette-wipe boundaries change. | Fractional speed, multiple within-frame emissions and trail replay improve motion representation; emission-time sampling by output column is absent. |
| [Fishbowl](../../effects/Fishbowl.h) | Anti-aliased warped multiline from an orientation trail. Head motion, noise deformation, trail phase/age and moving palette coordinates change. | Deep tweening and adaptive geometric midpoints refine the trail spatially; its noise and rendering state remain frame-wide. |
| [HopfFibration](../../effects/HopfFibration.h) | Projected fiber trails rasterized as polylines. S3 tumble phases, projected live tips, orientation and trail ages change. | Unoriented vector trails retain past fiber positions and are oriented at render time; that render state is not sampled per column. |
| [RingSpin](../../effects/RingSpin.h) | Fused ring-group scans of orientation trails. Ring orientations, trail sampling/weights and thickness at the head/tail change. | Orientation history, deep tweening and motion-blur rings are already present. Optional driver blanking changes emission duty, not animation sampling time. |

## Uncovered: particles

| Effect | Current rendering and changing state | Existing history or related machinery |
|---|---|---|
| [Galaxies](../../effects/Galaxies.h) | Forward point particles around six galaxy cores. Particle positions, velocities, spawn/kill, age opacity and frame shading change after one physics step. | Point-particle mode has no trail history; the pool capacity is 12,000, so per-particle state retention has material memory cost. |
| [MindSplatter](../../effects/MindSplatter.h) | Forward particle trails through fused or ordinary plotting paths. Physics positions, live tips, spawn/kill, trail ages, Mobius warp, hole gates and palette state change. | Stored trail anchors and an appended live tip support continuous-looking strokes; neither the warp nor the stroke state uses destination-column time. |

## Uncovered: feedback and lattice simulation

| Effect | Current rendering and changing state | Existing history or related machinery |
|---|---|---|
| [MeshFeedback](../../effects/MeshFeedback.h) | Oriented wireframe injected into a pixel feedback filter. Source rotation, mesh selection/shade and feedback transport, fading, noise/phase controls change. | begin_frame updates persistent feedback before plotting the source mesh. Feedback persistence is not a pair of temporal endpoints; both source and feedback presentation are uncovered. |
| [BZReactionDiffusion](../../effects/BZReactionDiffusion.h) | Frame-wide lattice shader after multiple reaction/diffusion substeps. Three species concentrations and lattice orientation change. | Q16 species arrays and spatial supersampling describe the latest field; physics substeps are not output-column times. |
| [GSReactionDiffusion](../../effects/GSReactionDiffusion.h) | Frame-wide lattice shader after Gray-Scott substeps. Concentrations, pigment, orientation, shimmer/noise, dissolve and reseeding change. | Physics-history scratch, hot-node culling and palette caches accelerate a frame; they do not retain all previous rendered field/material state for column-time sampling. |

## Uncovered: composed exception and custom inverse renderers

| Effect | Current rendering and changing state | Existing history or related machinery |
|---|---|---|
| [MobiusGrid](../../effects/MobiusGrid.h) | Composed pullback with temporal sampling disabled. Camera/source/warp clocks, palette and animated Mobius lens coefficients change. | Ordinary lens interpolation holds the previous coefficients until the endpoint; the effect has not opted into the composed temporal path. |
| [SphericalHarmonics](../../effects/SphericalHarmonics.h) | Full-screen inverse harmonic shader. Orientation, harmonic pair, morph amount, amplitude and palette state are prepared once per draw. | Harmonic morphing is continuous over frames, but the prepared frame is shared by all columns. |
| [HyperLattice](../../effects/HyperLattice.h) | Full-screen inverse ray tracing. Camera flight, 3D/4D rotation phases, pattern/shell controls, gain and iridescent palette state change. | Specialized trace paths and crossing scratch operate on one prepared frame context. |
| [Voronoi](../../effects/Voronoi.h) | Full-screen nearest-site shader. Site positions advance before KD-tree construction and block candidate selection. Site-count changes reseed the scene. | Tree and coherence blocks refer to one set of site positions for the entire frame. |
| [Raymarch](../../effects/Raymarch.h) | Object-bounded inverse volume scans, submitted in a forward object loop. Camera, per-object rotations, spin, shape controls and palette/noise phases change. | Ray marching samples space within each volume; object transforms, shading resources and cull bounds are held during that volume's draw. |

## Covered: composed column-time renderers

Every effect below declares `TEMPORAL_SAMPLING = true` in its linked source.
The common [implementation contract](pullback_temporal_sampling.md) covers
orientation interpolation, resource clocks, continuous controls and retained
color resources. Physical firmware enables it; the simulator selects it with
POV fidelity. Fidelity off deliberately presents the current instantaneous frame.

| Effect | Rendering content covered by the common temporal path |
|---|---|
| [AlienBrain](../../effects/AlienBrain.h) | Glitch-folded grid with animated wave shear. |
| [AlienCore](../../effects/AlienCore.h) | Mirrored grid through a glitch lens. |
| [AlienOcean](../../effects/AlienOcean.h) | Folded grid with mirrored drift. |
| [AshCloud](../../effects/AshCloud.h) | Curl-displaced lattice across dodecahedral facets. |
| [ChromaticLichen](../../effects/ChromaticLichen.h) | Glitch-folded grid with curl displacement. |
| [CosmicEyeball](../../effects/CosmicEyeball.h) | Mirrored grid with path-length hue. |
| [GridSpace](../../effects/GridSpace.h) | Affine lattice contours. |
| [JewelMelt](../../effects/JewelMelt.h) | Curl-displaced octahedral lattice. |
| [KaleidoscopeFlowers](../../effects/KaleidoscopeFlowers.h) | Dodecahedral folded grids. |
| [KaleidoscopeHexBright](../../effects/KaleidoscopeHexBright.h) | Hexagonal-prism twin wave. |
| [KaleidoscopeHexOil](../../effects/KaleidoscopeHexOil.h) | Displaced hexagonal-prism spiral with path-length hue. |
| [KaleidoscopeHexSoft](../../effects/KaleidoscopeHexSoft.h) | Drifting folded twin wave. |
| [KaleidoscopeMandala](../../effects/KaleidoscopeMandala.h) | Wave-sheared dodecahedral grid. |
| [KaleidoscopeMelt](../../effects/KaleidoscopeMelt.h) | Curl-displaced dodecahedral lattice. |
| [KaleidoscopePentBright](../../effects/KaleidoscopePentBright.h) | Pentagonal-prism polar lattice. |
| [KaleidoscopeSmooth](../../effects/KaleidoscopeSmooth.h) | Mirrored dodecahedral grids. |
| [KaleidoscopeStainedGlass](../../effects/KaleidoscopeStainedGlass.h) | Vector-noise-warped dodecahedral grid. |
| [MermaidSkin](../../effects/MermaidSkin.h) | Curl-displaced folded grid. |

Discrete parameter changes, topology changes, projection singularities and
spatial aliasing need separate diagnosis even for these covered effects. This
inventory does not treat an opt-in as a universal visual certification.

## Optional interpreter path

| Effect | Coverage and source evidence |
|---|---|
| [ShaderChain](../../workbench/shader/chain_host.h) | Uncovered. When HS_ENABLE_CHAIN_INTERPRETER enables its roster entry, draw_frame advances clocks, creates one frame context, prepares one interpreted program and runs a frame-wide shader. Its animated operation parameters/clocks and color resources do not use the composed effects' temporal opt-in. Authored effects and their interpreted workbench equivalents therefore need separate coverage accounting. |

## Audit boundary

This inventory covers every authored HS_EFFECT_LIST entry and its conditional
ShaderChain extension. The counts do not include test-only effects, profile
harnesses, primitive demos, or arbitrary future workbench programs as additional
authored effects.

Simulator effect-selection history mixing was fixed separately by atomic
replacement of complete frames and response history. It does not change the
25 uncovered authored renderers above. A screenshot taken while paused also
cannot by itself distinguish spatial and temporal artifacts.

No new solutions, priorities or performance promises are selected by this audit.
The earlier [design draft](remaining_temporal_seams.md) is provisional; use this
inventory as the scope for the next design discussion.
