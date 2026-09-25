# Finding 2: IslamicStars visual comparison

Open [the standalone gallery](index.html) or [the three-frame contact sheet](comparison.png).

These are deterministic full-canvas WASM renderer captures, not photographs of the physical display. Both builds use the same seed, frame sequence, 288 x 144 resolution, and default transition speed 1. The measured hardware profile used transition speed 4; this visual comparison uses the default choreography.

Baseline `b4a9848ff` differs from candidate `fa92838f9` only by reverting finding 2 (face bounds and its regression test). The full cycle covered 7104 frames and all 23 shapes, ending at the first shape's next spawn. Image values convert linear RGB16 to standard sRGB8 without exposure adjustment. Changes below 8-bit display quantization may be absent from the visible-difference count.

Selected images include the highest total absolute sRGB difference for every shape and the eight highest-difference frames overall, deduplicated. This intentionally emphasizes the strongest differences, not typical frames. All per-frame measurements are in `all-frame-metrics.json`; selected-frame metrics are in `manifest.json`. Difference views multiply absolute sRGB channel differences by 16 and clamp at 255. They do not depict normal visual intensity. Enlargements use nearest-pixel sampling.

3658 of 7104 frames have at least one changed sRGB pixel. Average changed pixels per frame: 45.59 of 41,472. Maximum changed pixels in any frame: 1414.

The gallery sphere previews project these same full-canvas images for orientation; they do not simulate LEDs, optical persistence, or the physical segmented driver. No disposition or landing of finding 2 is implied.

## Rendering limitation found during visual review

The user observed more dangling edges with the candidate and no corresponding artifacts in the simulator. The gallery sphere is a filled nearest-pixel texture projection, not the simulator's antialiased round LED geometry. Its latitude and longitude mapping also differ from Daydream. The flat before/after pixels remain useful as a raw renderer comparison, but neither gallery view establishes how the change looks in the simulator or on hardware.

Finding 2 widens the scan bounds while preserving the approximate half-plane distance outside convex corners. This can reveal miter extensions previously clipped by the old bounds. Its regression test establishes coverage of that approximate distance field, not improved geometric silhouette. Finding 2 remains deferred; a future visual evaluation must use matching simulator rendering and compare against exact corner distance.
