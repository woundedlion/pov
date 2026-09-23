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
 * @brief Slot-menu display labels, warp parameter labels, and stable export
 *        spellings for every enumerated field.
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
    "math::NoiseBasis::SIMPLEX", "math::NoiseBasis::FBM3",
    "math::NoiseBasis::RIDGED3"};
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
/** @brief Index of a warp parameter name in the per-position name tables. */
enum WarpParamName : uint8_t {
  WARP_NAME_TRANSLATION_X,
  WARP_NAME_TRANSLATION_Y,
  WARP_NAME_ROTATION,
  WARP_NAME_SCALE_X,
  WARP_NAME_SCALE_Y,
  WARP_NAME_SHEAR,
  WARP_NAME_FREQUENCY,
  WARP_NAME_FIELD_ANGLE,
  WARP_NAME_CENTER_X,
  WARP_NAME_CENTER_Y,
  WARP_NAME_RADIUS,
  WARP_NAME_TURNS,
  WARP_NAME_VECTOR_ANGLE,
  WARP_NAME_CELL_X,
  WARP_NAME_CELL_Y,
  WARP_NAME_OFFSET_X,
  WARP_NAME_OFFSET_Y,
  WARP_NAME_RADIAL_SCALE,
  WARP_NAME_RADIAL_PHASE,
  WARP_NAME_ANGULAR_PHASE,
  WARP_NAME_EDGE_WIDTH,
  WARP_NAME_CENTER_ORBIT,
  WARP_NAME_COUNT,
};

inline constexpr const char *OUTER_WARP_PARAM_NAMES[] = {
    "Planar Warp 1 Translation X", "Planar Warp 1 Translation Y",
    "Planar Warp 1 Rotation",      "Planar Warp 1 Scale X",
    "Planar Warp 1 Scale Y",       "Planar Warp 1 Shear",
    "Planar Warp 1 Frequency",     "Planar Warp 1 Field Angle",
    "Planar Warp 1 Center X",      "Planar Warp 1 Center Y",
    "Planar Warp 1 Radius",        "Planar Warp 1 Turns",
    "Planar Warp 1 Vector Angle",  "Planar Warp 1 Cell X",
    "Planar Warp 1 Cell Y",        "Planar Warp 1 Offset X",
    "Planar Warp 1 Offset Y",      "Planar Warp 1 Radial Scale",
    "Planar Warp 1 Radial Phase",  "Planar Warp 1 Angular Phase",
    "Planar Warp 1 Edge Width",    "Planar Warp 1 Center Orbit"};
inline constexpr const char *INNER_WARP_PARAM_NAMES[] = {
    "Planar Warp 2 Translation X", "Planar Warp 2 Translation Y",
    "Planar Warp 2 Rotation",      "Planar Warp 2 Scale X",
    "Planar Warp 2 Scale Y",       "Planar Warp 2 Shear",
    "Planar Warp 2 Frequency",     "Planar Warp 2 Field Angle",
    "Planar Warp 2 Center X",      "Planar Warp 2 Center Y",
    "Planar Warp 2 Radius",        "Planar Warp 2 Turns",
    "Planar Warp 2 Vector Angle",  "Planar Warp 2 Cell X",
    "Planar Warp 2 Cell Y",        "Planar Warp 2 Offset X",
    "Planar Warp 2 Offset Y",      "Planar Warp 2 Radial Scale",
    "Planar Warp 2 Radial Phase",  "Planar Warp 2 Angular Phase",
    "Planar Warp 2 Edge Width",    "Planar Warp 2 Center Orbit"};
static_assert(sizeof(OUTER_WARP_PARAM_NAMES) / sizeof(const char *) ==
                  WARP_NAME_COUNT,
              "outer warp name table must match WarpParamName");
static_assert(sizeof(INNER_WARP_PARAM_NAMES) / sizeof(const char *) ==
                  WARP_NAME_COUNT,
              "inner warp name table must match WarpParamName");
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

static_assert(std::size(FUNCTION_OPTIONS) ==
              std::size(FUNCTION_EXPORT_OPTIONS));
static_assert(std::size(TESSELLATION_KIND_OPTIONS) ==
              std::size(TESSELLATION_KIND_EXPORT_OPTIONS));
static_assert(std::size(PROJECTION_OPTIONS) ==
              std::size(PROJECTION_EXPORT_OPTIONS));
static_assert(std::size(PEIRCE_LAYOUT_OPTIONS) ==
              std::size(PEIRCE_LAYOUT_EXPORT_OPTIONS));
static_assert(std::size(AIROCEAN_LAYOUT_OPTIONS) ==
              std::size(AIROCEAN_LAYOUT_EXPORT_OPTIONS));
static_assert(std::size(BONNE_HEMISPHERE_OPTIONS) ==
              std::size(BONNE_HEMISPHERE_EXPORT_OPTIONS));
static_assert(std::size(GNOMONIC_HEMISPHERE_OPTIONS) ==
              std::size(GNOMONIC_HEMISPHERE_EXPORT_OPTIONS));
static_assert(std::size(PROJECTION_FRAME_OPTIONS) ==
              std::size(PROJECTION_FRAME_EXPORT_OPTIONS));
static_assert(std::size(LENS_OPTIONS) == std::size(LENS_EXPORT_OPTIONS));
static_assert(std::size(SURFACE_NOISE_OPTIONS) ==
              std::size(SURFACE_NOISE_EXPORT_OPTIONS));
static_assert(std::size(SURFACE_NOISE_PLACEMENT_OPTIONS) ==
              std::size(SURFACE_NOISE_PLACEMENT_EXPORT_OPTIONS));
static_assert(std::size(SURFACE_CURL_INTEGRATOR_OPTIONS) ==
              std::size(SURFACE_CURL_INTEGRATOR_EXPORT_OPTIONS));
static_assert(std::size(WARP_OPTIONS) == std::size(WARP_EXPORT_OPTIONS));
static_assert(std::size(NOISE_BASIS_OPTIONS) ==
              std::size(NOISE_BASIS_EXPORT_OPTIONS));
static_assert(std::size(POLAR_MODE_OPTIONS) ==
              std::size(POLAR_MODE_EXPORT_OPTIONS));
static_assert(std::size(CURL_INTEGRATOR_OPTIONS) ==
              std::size(CURL_INTEGRATOR_EXPORT_OPTIONS));
static_assert(std::size(WARP_ENVELOPE_OPTIONS) ==
              std::size(WARP_ENVELOPE_EXPORT_OPTIONS));
static_assert(std::size(SIGNAL_OPTIONS) == std::size(SIGNAL_EXPORT_OPTIONS));
static_assert(std::size(VALUE_TRANSFER_OPTIONS) ==
              std::size(VALUE_TRANSFER_EXPORT_OPTIONS));
static_assert(std::size(COVERAGE_OPTIONS) ==
              std::size(COVERAGE_EXPORT_OPTIONS));
static_assert(std::size(PALETTE_OPTIONS) == std::size(PALETTE_EXPORT_OPTIONS));
static_assert(std::size(PALETTE_MAPPING_OPTIONS) ==
              std::size(PALETTE_MAPPING_EXPORT_OPTIONS));
static_assert(std::size(BRIGHTNESS_ENVELOPE_OPTIONS) ==
              std::size(BRIGHTNESS_ENVELOPE_EXPORT_OPTIONS));
static_assert(std::size(HUE_SHIFT_OPTIONS) ==
              std::size(HUE_SHIFT_EXPORT_OPTIONS));

/** @brief Length of a table whose last row is @p last, the enum's final
    in-range value. */
template <typename Enum> inline constexpr int option_count_through(Enum last) {
  return static_cast<int>(last) + 1;
}
static_assert(NUM_FUNCTIONS == option_count_through(Function::TESSELLATION));
static_assert(
    NUM_TESSELLATION_KINDS ==
    option_count_through(Pullback::Source::TessellationKind::HEXAGONAL));
static_assert(NUM_PROJECTIONS ==
              option_count_through(Projection::EQUIRECTANGULAR));
static_assert(NUM_PEIRCE_LAYOUTS ==
              option_count_through(PeirceLayout::VERTICAL));
static_assert(NUM_AIROCEAN_LAYOUTS ==
              option_count_through(AiroceanLayout::HORIZONTAL));
static_assert(NUM_BONNE_HEMISPHERES ==
              option_count_through(BonneHemisphere::SOUTH));
static_assert(NUM_GNOMONIC_HEMISPHERES ==
              option_count_through(GnomonicHemispherePolicy::BACK_HEMISPHERE));
static_assert(NUM_PROJECTION_FRAMES ==
              option_count_through(ProjectionFramePolicy::SPIN_WANDER));
static_assert(NUM_LENSES ==
              option_count_through(SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM));
static_assert(NUM_SURFACE_NOISE == option_count_through(SurfaceNoise::CURL));
static_assert(NUM_SURFACE_NOISE_PLACEMENTS ==
              option_count_through(SurfaceNoisePlacement::AFTER_LENS));
static_assert(NUM_SURFACE_CURL_INTEGRATORS ==
              option_count_through(SurfaceCurlIntegrator::MIDPOINT_2X));
static_assert(NUM_WARPS == option_count_through(WarpStageKind::POLAR_CHART));
static_assert(NUM_NOISE_BASES ==
              option_count_through(math::NoiseBasis::RIDGED3));
static_assert(NUM_POLAR_MODES == option_count_through(PolarMode::LOGARITHMIC));
static_assert(NUM_CURL_INTEGRATORS ==
              option_count_through(CurlIntegrator::MIDPOINT_4));
static_assert(NUM_WARP_ENVELOPES ==
              option_count_through(WarpEnvelope::EDGE_FADE));
static_assert(NUM_SIGNALS == option_count_through(SignalWeight::PROJECTION));
static_assert(NUM_VALUE_TRANSFERS ==
              option_count_through(ValueTransfer::SMOOTH_BANDS));
static_assert(NUM_COVERAGE_POLICIES ==
              option_count_through(CoveragePolicy::PROJECTION_WEIGHT));
static_assert(NUM_PALETTES == option_count_through(PaletteMode::ANALOGOUS));
static_assert(NUM_PALETTE_MAPPINGS ==
              option_count_through(PaletteMapping::REVERSE));
static_assert(NUM_BRIGHTNESS_ENVELOPES ==
              option_count_through(BrightnessEnvelope::DESCENDING));
static_assert(NUM_HUE_SHIFT_MODES ==
              option_count_through(HueShiftMode::WARP_DISPLACEMENT));

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
