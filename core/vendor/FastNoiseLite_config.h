// FastNoiseLite configuration for Holosphere.
// This macro hard-routes GenNoiseSingle straight to OpenSimplex2 (bypassing the
// per-noise-type switch) and skips the 3D warp-transform setup. It does not
// itself remove the Cellular/Perlin/Value/OpenSimplex2S, fractal, or domain-warp
// code — those definitions stay in the header and are dropped by the compiler as
// dead code once nothing references them, which is what yields the lean binary.
// Remove this define to re-enable the full FastNoiseLite feature set.
#define FASTNOISELITE_ONLY_OPENSIMPLEX2
