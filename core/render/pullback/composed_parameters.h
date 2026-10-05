/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/**
 * @brief Interpolates one parameter family across a preset transition.
 * @details Driven by the family's field table; each field moves on the curve
 * its descriptor names, and every member the table does not cover, the Mobius
 * coefficients included, snaps to @p b once progress reaches 1.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param t Progress fraction.
 * @return The interpolated family.
 */
template <Pullback::HasFields T>
inline T interpolate(const T &a, const T &b, float t) {
  return Pullback::Fields::interpolate(a, b, t);
}

inline MobiusLensParams interpolate(const MobiusLensParams &a,
                                    const MobiusLensParams &b, float t) {
  MobiusLensParams value;
  value.mobius = t < 1.0f ? a.mobius : b.mobius;
  return value;
}

/**
 * @brief Interpolates a whole parameter set, family by family.
 * @param from Parameters at progress 0.
 * @param to Parameters at progress 1.
 * @param progress Progress fraction, already eased by the caller.
 * @return The interpolated parameter set.
 */
template <typename... Resources>
inline ComposedDetail::ParameterSet<ComposedDetail::ResourceList<Resources...>>
interpolate(const ComposedDetail::ParameterSet<
                ComposedDetail::ResourceList<Resources...>> &from,
            const ComposedDetail::ParameterSet<
                ComposedDetail::ResourceList<Resources...>> &to,
            float progress) {
  return {ComposedDetail::ParameterBlock<Resources::KEY,
                                         typename Resources::Family>{
      interpolate(from.template get<Resources::KEY>(),
                  to.template get<Resources::KEY>(), progress)}...};
}

/**
 * @brief Whether every field of a parameter family is inside its authored
 *        range.
 * @details Driven by the family's field table, so the admissibility ranges a
 * restored snapshot must pass are the same descriptors the sliders register
 * with.
 */
template <Pullback::HasFields T> inline bool valid(const T &value) {
  return Pullback::Fields::valid(value);
}

inline bool valid(const MobiusLensParams &p) {
  const float values[] = {p.mobius.a.re, p.mobius.a.im, p.mobius.b.re,
                          p.mobius.b.im, p.mobius.c.re, p.mobius.c.im,
                          p.mobius.d.re, p.mobius.d.im};
  for (float value : values)
    if (!std::isfinite(value) ||
        fabsf(value) > MobiusLensParams::COEFFICIENT_LIMIT)
      return false;
  return MobiusLensParams::nondegenerate(p.mobius);
}

inline bool valid(const ColorParams &p) {
  return Pullback::Fields::valid(p) &&
         static_cast<uint8_t>(p.palette_mapping) <=
             static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE);
}

/**
 * @brief Whether every family of a parameter set is in range.
 * @return True only when every declared instance passes.
 */
template <typename... Resources>
inline bool valid(const ComposedDetail::ParameterSet<
                  ComposedDetail::ResourceList<Resources...>> &params) {
  bool result = true;
  params.visit([&]<typename Resource>(const auto &family) {
    result = valid(family) && result;
  });
  return result;
}

template <bool Enabled> struct OptionalNoise {};
template <> struct OptionalNoise<true> {
  FastNoiseLite noise;
};

/** @brief Hue-rotation LUT storage; empty when the effect never rotates hue. */
template <bool Enabled> struct OptionalHueRotationLut {};
template <> struct OptionalHueRotationLut<true> {
  std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE> hue_rotation_lut;
  /** Palette bake the resident table was built from; 0 matches no bake,
      forcing the first build. */
  uint32_t hue_rotation_lut_bake = 0;
};

/** @brief Hue-noise LUT and the inputs it was baked from; empty unless the
    hue source is the noise field. */
template <bool Enabled> struct OptionalHueNoiseLut {};
template <> struct OptionalHueNoiseLut<true> {
  FastNoiseLite color_noise;
  std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut;
  Pullback::Color::HueNoiseBakeCache hue_noise_bake;
};

/** @brief Projection-walk noise storage; empty when disabled. */
template <bool Enabled> struct ProjectionWalkNoise {};
template <> struct ProjectionWalkNoise<true> {
  FastNoiseLite projection_walk_noise;
};

/** @brief Persistent projection-walk state; empty when disabled. */
template <bool Enabled> struct ProjectionWalkState {
  math::Quaternion frame_conjugate() const { return math::Quaternion(); }
};
template <> struct ProjectionWalkState<true> {
  math::Orientation<> projection_walk;
  math::Quaternion projection_walk_previous;
  math::Quaternion projection_wander;
  math::Quaternion projection_conjugate;
  math::Quaternion base_orientation = Pullback::projection_base_orientation();
  float projection_spin = 0.0f;

  math::Quaternion frame_conjugate() const { return projection_conjugate; }
};

/** @brief Sphere-to-plane projection of a composed effect's Stage::Project. */
enum class ProjectionKind : uint8_t {
  STEREOGRAPHIC,
  GNOMONIC_FOLDED,
  EQUIRECTANGULAR,
  FOLDED_SINUSOIDAL
};

/** @brief Whether @p projection reads the central-meridian field. */
constexpr bool uses_central_meridian(ProjectionKind projection) {
  return projection == ProjectionKind::EQUIRECTANGULAR ||
         projection == ProjectionKind::FOLDED_SINUSOIDAL;
}

/** @brief Whether @p projection reads the singularity-fade field. Folded
    sinusoidal has no singular locus and returns fixed weights. */
constexpr bool uses_singularity_fade(ProjectionKind projection) {
  return projection != ProjectionKind::FOLDED_SINUSOIDAL;
}

/** @brief Optional transfer curve an effect's material stage composes. */
enum class TransferKind : uint8_t { NONE, ISO_CONTOUR };

/** @brief Optional value-dependent coverage stage after sampling. */
enum class FieldCoverageKind : uint8_t { NONE, VALUE_CUTOUT };
