# Finding 21 visual comparison

These are illustrative native CPU renders at the production canvas resolution,
288 by 144. They are not photographs, Teensy framebuffer captures, or proof of
visibility on the rotating display.

Candidate source: `d4b8d7b75f3e2a87ccd2801a946d0b8ad57e6cc6`.
Baseline source behavior: `0c02f3912`, reproduced by disabling only the new
vertex-cone block in `Face::plane_dist_convex`. The original half-plane loop and
all other rendering code run unchanged in the same diagnostic executable.
The temporary harness patch, source and raw RGB16 captures are retained locally
in `C:/work/Holosphere/.git/implement-sessions/face21-visuals/`. No diagnostic edits
are committed to production source.

Each before/after effect pass resets the engine's global state, arenas,
timeline, display buffers and RNG (seed 1337). The deterministic mock clock
advances 33 milliseconds per frame. Frame numbers here are zero-based.
Rendering covers the full 288x144 canvas. Device profiling uses a 62.5 ms
display window and a segmented viewport, so similarly numbered host frames
are representative illustrations, not exact copies of captured device poses.

The default-schedule sweep captures frames 0 through 720 inclusive at intervals
of 12: 61 samples for HankinSolids and 61 for IslamicStars. The extended sweep
uses `HS_PROFILE_ORDERED_CYCLE=1` and `HS_PROFILE_TRANS_SPEED=4`, captures every
120 frames, plus the neighborhood of zero-based Hankin frame 2251 and Islamic
frame 2810. It ends at Hankin frame 2280 and Islamic frame 2840.

Raw `.rgb16` files contain row-major little-endian unsigned 16-bit RGB channels
in linear light. PNG conversion uses the repository's complete
`linear_to_srgb_lut`, whose outputs are identical to `linear_to_srgb8`. The same
fixed mapping is used for every image: no exposure adjustment or normalization.
Native comparison panels display the 288x144 images at 1:1. Zoom panels use
nearest-neighbor spatial enlargement by eight. Difference panels show absolute
8-bit sRGB differences with an explicitly labeled gain of eight, clipped at
255; standalone `diff_1x` PNGs have no gain.

The synthetic test is a white triangle with tangent-plane vertices
(0.20,0), (-0.10,0.015), (-0.10,-0.015), about a 5.7-degree apex angle. Its axis
has azimuth 1.31 radians and polar angle `0.38 + config * 0.12` radians for
configurations 0 through 19. These deliberately acute cases expose the
half-plane approximation's elongated tip. They are not claimed to represent
the frequency or contrast of corners in the effects.

`metrics.json` contains every captured matching pair and complete per-group
distributions. Changed-pixel counts use the maximum absolute channel change
per pixel. Thresholds 1, 8 and 16 are on the 0..255 sRGB scale. Linear-light
metrics use the original 0..65535 channels. `native_comparison` and
`zoom_comparison` select the sample with the greatest total absolute linear
channel difference. `peakpixel_comparison` separately selects the sample with
the greatest individual sRGB channel difference; it can differ substantially
from the sample with the largest total change. These selections are visual
stress cases among sampled frames, not statistical estimates of all frames.

The candidate leaves cull bounds unchanged and raises the exterior distance
near convex vertices. It removes inflated corner coverage; it does not recover
light at previously culled pixels. The white-triangle comparison dims pixels
only. Color differences in complete effects may brighten some channels due to
changed overlap or distance-dependent shading.

At native resolution the tested effect comparisons are difficult to distinguish
side by side. Differences become clearer in amplified difference panels and
enlarged corner crops. Some isolated pixels have substantial contrast changes;
the evidence does not support claiming that every difference is invisible.
The synthetic acute white triangle has a clearly shortened tip. These host
images do not demonstrate a large overall visual improvement in the tested
effects, or establish visibility on the physical display.
