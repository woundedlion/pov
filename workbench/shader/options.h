/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file options.h
 * @brief Slot-menu vocabulary: the display label and stable export spelling of every enumerated field.
 */

#include "workbench/shader/config.h"

namespace Workbench {

inline constexpr const char *FUNCTION_OPTIONS[] = {"Twin Wave",
                                                   "Rings",
                                                   "Spiral",
                                                   "Grid",
                                                   "Noise Contour (Projected)",
                                                   "Primitive Lattice",
                                                   "Noise Contour (Sphere)",
                                                   "Spherical Rings",
                                                   "Escape Fractal",
                                                   "Tessellation"};
inline constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
    "Function::TWIN_WAVE",
    "Function::RINGS",
    "Function::SPIRAL",
    "Function::GRID",
    "Function::NOISE_CONTOUR",
    "Function::PRIMITIVE_LATTICE",
    "Function::NOISE_CONTOUR_SPHERE",
    "Function::SPHERICAL_RINGS",
    "Function::FRACTAL",
    "Function::TESSELLATION"};
inline constexpr int NUM_FUNCTIONS = std::size(FUNCTION_OPTIONS);
inline constexpr const char *TESSELLATION_KIND_OPTIONS[] = {
    "Triangular", "Square", "Hexagonal"};
inline constexpr const char *TESSELLATION_KIND_EXPORT_OPTIONS[] = {
    "Pullback::Source::TessellationKind::TRIANGULAR",
    "Pullback::Source::TessellationKind::SQUARE",
    "Pullback::Source::TessellationKind::HEXAGONAL"};
inline constexpr int NUM_TESSELLATION_KINDS =
    std::size(TESSELLATION_KIND_OPTIONS);
inline constexpr const char *PROJECTION_OPTIONS[] = {
    "Folded Sinusoidal",  "Stereographic",       "Gnomonic",       "Bonne",
    "Peirce Quincuncial", "Dymaxion / Airocean", "Equirectangular"};
inline constexpr const char *PROJECTION_EXPORT_OPTIONS[] = {
    "Projection::SINUSOIDAL",         "Projection::STEREOGRAPHIC",
    "Projection::GNOMONIC",           "Projection::BONNE",
    "Projection::PEIRCE_QUINCUNCIAL", "Projection::AIROCEAN",
    "Projection::EQUIRECTANGULAR"};
inline constexpr int NUM_PROJECTIONS = std::size(PROJECTION_OPTIONS);
inline constexpr const char *PEIRCE_LAYOUT_OPTIONS[] = {
    "Diamond", "Square", "Horizontal", "Vertical"};
inline constexpr const char *PEIRCE_LAYOUT_EXPORT_OPTIONS[] = {
    "PeirceLayout::DIAMOND", "PeirceLayout::SQUARE", "PeirceLayout::HORIZONTAL",
    "PeirceLayout::VERTICAL"};
inline constexpr int NUM_PEIRCE_LAYOUTS = std::size(PEIRCE_LAYOUT_OPTIONS);
inline constexpr const char *AIROCEAN_LAYOUT_OPTIONS[] = {"Vertical",
                                                          "Horizontal"};
inline constexpr const char *AIROCEAN_LAYOUT_EXPORT_OPTIONS[] = {
    "AiroceanLayout::VERTICAL", "AiroceanLayout::HORIZONTAL"};
inline constexpr int NUM_AIROCEAN_LAYOUTS = std::size(AIROCEAN_LAYOUT_OPTIONS);
inline constexpr const char *BONNE_HEMISPHERE_OPTIONS[] = {"North", "South"};
inline constexpr const char *BONNE_HEMISPHERE_EXPORT_OPTIONS[] = {
    "BonneHemisphere::NORTH", "BonneHemisphere::SOUTH"};
inline constexpr int NUM_BONNE_HEMISPHERES =
    std::size(BONNE_HEMISPHERE_OPTIONS);
inline constexpr const char *GNOMONIC_HEMISPHERE_OPTIONS[] = {
    "Folded", "Front Hemisphere", "Back Hemisphere"};
inline constexpr const char *GNOMONIC_HEMISPHERE_EXPORT_OPTIONS[] = {
    "GnomonicHemispherePolicy::FOLDED",
    "GnomonicHemispherePolicy::FRONT_HEMISPHERE",
    "GnomonicHemispherePolicy::BACK_HEMISPHERE"};
inline constexpr int NUM_GNOMONIC_HEMISPHERES =
    std::size(GNOMONIC_HEMISPHERE_OPTIONS);
inline constexpr const char *PROJECTION_FRAME_OPTIONS[] = {"Identity",
                                                           "Spin + Wander"};
inline constexpr const char *PROJECTION_FRAME_EXPORT_OPTIONS[] = {
    "ProjectionFramePolicy::IDENTITY", "ProjectionFramePolicy::SPIN_WANDER"};
inline constexpr int NUM_PROJECTION_FRAMES =
    std::size(PROJECTION_FRAME_OPTIONS);
inline constexpr const char *LENS_OPTIONS[] = {
    "None",
    "Glitch",
    "Twist",
    "Kaleidoscope (Azimuthal 6-fold)",
    "Mobius",
    "Kaleidoscope (Tetrahedral)",
    "Kaleidoscope (Octahedral / Cubic)",
    "Kaleidoscope (Dodecahedral / Icosahedral)",
    "Kaleidoscope (Triangular Prism)",
    "Kaleidoscope (Square Prism)",
    "Kaleidoscope (Pentagonal Prism)",
    "Kaleidoscope (Hexagonal Prism)",
    "Kaleidoscope (Octagonal Prism)"};
inline constexpr const char *LENS_EXPORT_OPTIONS[] = {
    "SurfaceLens::NONE",
    "SurfaceLens::GLITCH",
    "SurfaceLens::TWIST",
    "SurfaceLens::KALEIDOSCOPE",
    "SurfaceLens::MOBIUS",
    "SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL",
    "SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL",
    "SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL",
    "SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM",
    "SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM",
    "SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM",
    "SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM",
    "SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM"};
inline constexpr int NUM_LENSES = std::size(LENS_OPTIONS);
inline constexpr const char *SURFACE_NOISE_OPTIONS[] = {"None", "Direct",
                                                        "Curl"};
inline constexpr const char *SURFACE_NOISE_EXPORT_OPTIONS[] = {
    "SurfaceNoise::NONE", "SurfaceNoise::DIRECT", "SurfaceNoise::CURL"};
inline constexpr int NUM_SURFACE_NOISE = std::size(SURFACE_NOISE_OPTIONS);
inline constexpr const char *SURFACE_NOISE_PLACEMENT_OPTIONS[] = {"Before Lens",
                                                                  "After Lens"};
inline constexpr const char *SURFACE_NOISE_PLACEMENT_EXPORT_OPTIONS[] = {
    "SurfaceNoisePlacement::BEFORE_LENS", "SurfaceNoisePlacement::AFTER_LENS"};
inline constexpr int NUM_SURFACE_NOISE_PLACEMENTS =
    std::size(SURFACE_NOISE_PLACEMENT_OPTIONS);
inline constexpr const char *SURFACE_CURL_INTEGRATOR_OPTIONS[] = {
    "Euler", "Midpoint", "Midpoint 2x"};
inline constexpr const char *SURFACE_CURL_INTEGRATOR_EXPORT_OPTIONS[] = {
    "SurfaceCurlIntegrator::EULER", "SurfaceCurlIntegrator::MIDPOINT",
    "SurfaceCurlIntegrator::MIDPOINT_2X"};
inline constexpr int NUM_SURFACE_CURL_INTEGRATORS =
    std::size(SURFACE_CURL_INTEGRATOR_OPTIONS);
inline constexpr const char *WARP_OPTIONS[] = {"None",
                                               "Affine Frame",
                                               "Wave Shear",
                                               "Vortex",
                                               "Projected Vector Noise",
                                               "Projected Curl Flow",
                                               "Mirror Tile",
                                               "Polar Chart"};
inline constexpr const char *WARP_EXPORT_OPTIONS[] = {
    "WarpStageKind::NONE",         "WarpStageKind::AFFINE_FRAME",
    "WarpStageKind::WAVE_SHEAR",   "WarpStageKind::VORTEX",
    "WarpStageKind::VECTOR_NOISE", "WarpStageKind::CURL_FLOW",
    "WarpStageKind::MIRROR_TILE",  "WarpStageKind::POLAR_CHART"};
inline constexpr int NUM_WARPS = std::size(WARP_OPTIONS);
inline constexpr const char *warp_option(WarpStageKind kind) {
  const uint8_t index = static_cast<uint8_t>(kind);
  return index < NUM_WARPS ? WARP_OPTIONS[index] : "Legacy Stereo Noise";
}
inline constexpr const char *NOISE_BASIS_OPTIONS[] = {"Simplex", "FBM 3",
                                                      "Ridged 3"};
inline constexpr const char *NOISE_BASIS_EXPORT_OPTIONS[] = {
    "NoiseBasis::SIMPLEX", "NoiseBasis::FBM3", "NoiseBasis::RIDGED3"};
inline constexpr int NUM_NOISE_BASES = std::size(NOISE_BASIS_OPTIONS);
inline constexpr const char *POLAR_MODE_OPTIONS[] = {"Linear", "Logarithmic"};
inline constexpr const char *POLAR_MODE_EXPORT_OPTIONS[] = {
    "PolarMode::LINEAR", "PolarMode::LOGARITHMIC"};
inline constexpr int NUM_POLAR_MODES = std::size(POLAR_MODE_OPTIONS);
inline constexpr const char *CURL_INTEGRATOR_OPTIONS[] = {
    "Euler 1", "Midpoint 2", "Midpoint 4"};
inline constexpr const char *CURL_INTEGRATOR_EXPORT_OPTIONS[] = {
    "CurlIntegrator::EULER_1", "CurlIntegrator::MIDPOINT_2",
    "CurlIntegrator::MIDPOINT_4"};
inline constexpr int NUM_CURL_INTEGRATORS = std::size(CURL_INTEGRATOR_OPTIONS);
inline constexpr const char *WARP_ENVELOPE_OPTIONS[] = {
    "Flat", "Projection Weight", "Edge Fade"};
inline constexpr const char *WARP_ENVELOPE_EXPORT_OPTIONS[] = {
    "WarpEnvelope::FLAT", "WarpEnvelope::PROJECTION_WEIGHT",
    "WarpEnvelope::EDGE_FADE"};
inline constexpr int NUM_WARP_ENVELOPES = std::size(WARP_ENVELOPE_OPTIONS);
inline constexpr const char *SIGNAL_OPTIONS[] = {"None", "Projection"};
inline constexpr const char *SIGNAL_EXPORT_OPTIONS[] = {
    "SignalWeight::NONE", "SignalWeight::PROJECTION"};
inline constexpr int NUM_SIGNALS = std::size(SIGNAL_OPTIONS);
inline constexpr const char *VALUE_TRANSFER_OPTIONS[] = {
    "None", "Ridge", "Iso Contour", "Smooth Bands"};
inline constexpr const char *VALUE_TRANSFER_EXPORT_OPTIONS[] = {
    "ValueTransfer::NONE", "ValueTransfer::RIDGE", "ValueTransfer::ISO_CONTOUR",
    "ValueTransfer::SMOOTH_BANDS"};
inline constexpr int NUM_VALUE_TRANSFERS = std::size(VALUE_TRANSFER_OPTIONS);
inline constexpr const char *COVERAGE_OPTIONS[] = {
    "Opaque", "Projection Weight Squared", "Value Cutout", "Edge Fade",
    "Projection Weight"};
inline constexpr const char *COVERAGE_EXPORT_OPTIONS[] = {
    "CoveragePolicy::OPAQUE", "CoveragePolicy::PROJECTION_WEIGHT_SQUARED",
    "CoveragePolicy::VALUE_CUTOUT", "CoveragePolicy::EDGE_FADE",
    "CoveragePolicy::PROJECTION_WEIGHT"};
inline constexpr int NUM_COVERAGE_POLICIES = std::size(COVERAGE_OPTIONS);
inline constexpr const char *PALETTE_OPTIONS[] = {
    "Generated Triadic", "Generated Complementary", "Generated Analogous"};
inline constexpr const char *PALETTE_EXPORT_OPTIONS[] = {
    "PaletteMode::TRIADIC", "PaletteMode::COMPLEMENTARY",
    "PaletteMode::ANALOGOUS"};
inline constexpr int NUM_PALETTES = std::size(PALETTE_OPTIONS);
inline constexpr const char *PALETTE_MAPPING_OPTIONS[] = {"Cup", "Bell",
                                                          "Linear", "Reverse"};
inline constexpr const char *PALETTE_MAPPING_EXPORT_OPTIONS[] = {
    "PaletteMapping::CUP", "PaletteMapping::BELL", "PaletteMapping::LINEAR",
    "PaletteMapping::REVERSE"};
inline constexpr int NUM_PALETTE_MAPPINGS = std::size(PALETTE_MAPPING_OPTIONS);
inline constexpr const char *BRIGHTNESS_ENVELOPE_OPTIONS[] = {
    "None", "Cup", "Bell", "Ascending", "Descending"};
inline constexpr const char *BRIGHTNESS_ENVELOPE_EXPORT_OPTIONS[] = {
    "BrightnessEnvelope::NONE", "BrightnessEnvelope::CUP",
    "BrightnessEnvelope::BELL", "BrightnessEnvelope::ASCENDING",
    "BrightnessEnvelope::DESCENDING"};
inline constexpr int NUM_BRIGHTNESS_ENVELOPES =
    std::size(BRIGHTNESS_ENVELOPE_OPTIONS);
inline constexpr const char *HUE_SHIFT_OPTIONS[] = {"None", "Noise",
                                                    "Total Warp Displacement"};
inline constexpr const char *HUE_SHIFT_EXPORT_OPTIONS[] = {
    "HueShiftMode::NONE", "HueShiftMode::NOISE",
    "HueShiftMode::WARP_DISPLACEMENT"};
inline constexpr int NUM_HUE_SHIFT_MODES = std::size(HUE_SHIFT_OPTIONS);

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
