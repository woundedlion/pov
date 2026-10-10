/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_parameters.h
 * @brief Parameter-family interpolation, range validation and optional
 * runtime storage for composed effects.
 */

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

/**
 * @brief Mobius lens transition: coefficients hold @p a, then snap to @p b.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param t Progress fraction.
 * @return @p a while @p t < 1, else @p b.
 */
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
 * @details Ranges come from the family's field table.
 * @tparam T Parameter family with a field table.
 * @param value Family to check.
 * @return True when every field is in range.
 */
template <Pullback::HasFields T> inline bool valid(const T &value) {
  return Pullback::Fields::valid(value);
}

/**
 * @brief Whether Mobius coefficients are finite, within
 *        `MobiusLensParams::COEFFICIENT_LIMIT` and nondegenerate.
 * @param p Lens family to check.
 * @return True when the lens is usable.
 */
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

/**
 * @brief Whether a colour family's fields and palette mapping are in range.
 * @param p Colour family to check.
 * @return True when every field and the mapping enumerator are valid.
 */
inline bool valid(const ColorParams &p) {
  return Pullback::Fields::valid(p) &&
         static_cast<uint8_t>(p.palette_mapping) <=
             static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE);
}

/**
 * @brief Whether every family of a parameter set is in range.
 * @tparam Resources The set's parameter resources.
 * @param params Parameter set to check.
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

/** @brief Noise-field storage; empty when disabled.
    @tparam Enabled Whether the field is stored. */
template <bool Enabled> struct OptionalNoise {};
/** @brief Enabled noise-field storage. */
template <> struct OptionalNoise<true> {
  FastNoiseLite noise; ///< The resource's noise field.
};

/** @brief Hue-rotation LUT storage; empty when the effect never rotates hue. */
template <bool Enabled> struct OptionalHueRotationLut {};
/** @brief Enabled hue-rotation LUT storage. */
template <> struct OptionalHueRotationLut<true> {
  /// Palette colours by palette coordinate and hue-rotation step.
  std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE> hue_rotation_lut;
  /** Palette bake the resident table was built from; 0 matches no bake,
      forcing the first build. */
  uint32_t hue_rotation_lut_bake = 0;
};

/** @brief Hue-noise LUT and the inputs it was baked from; empty unless the
    hue source is the noise field. */
template <bool Enabled> struct OptionalHueNoiseLut {};
/** @brief Enabled hue-noise storage. */
template <> struct OptionalHueNoiseLut<true> {
  FastNoiseLite color_noise; ///< Hue-noise field the LUT is baked from.
  /// Cube-face hue-noise samples.
  std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut;
  /// Scale and phase the resident LUT was baked from.
  Pullback::Color::HueNoiseBakeCache hue_noise_bake;
};

/** @brief Projection-walk noise storage; empty when disabled. */
template <bool Enabled> struct ProjectionWalkNoise {};
/** @brief Enabled projection-walk noise storage. */
template <> struct ProjectionWalkNoise<true> {
  FastNoiseLite projection_walk_noise; ///< Noise driving the projection walk.
};

/** @brief Persistent projection-walk state; empty when disabled. */
template <bool Enabled> struct ProjectionWalkState {
  /** @brief Inverse projection frame rotation.
      @return The identity: a static projection never rotates. */
  math::Quaternion frame_conjugate() const { return math::Quaternion(); }
};
/** @brief Enabled projection-walk state. */
template <> struct ProjectionWalkState<true> {
  math::Orientation<> projection_walk; ///< Random-walk orientation.
  /** Walk orientation at the previous frame, for the per-frame delta. */
  math::Quaternion projection_walk_previous;
  /** Accumulated walk deltas, each scaled by the wander parameter. */
  math::Quaternion projection_wander;
  /** Inverse of spin * `base_orientation` * `projection_wander`. */
  math::Quaternion projection_conjugate;
  /** Fixed orientation the spin and wander compose onto. */
  math::Quaternion base_orientation = Pullback::projection_base_orientation();
  float projection_spin = 0.0f; ///< Spin about Y, radians, |value| < 2 pi.

  /** @brief Inverse projection frame rotation.
      @return `projection_conjugate`. */
  math::Quaternion frame_conjugate() const { return projection_conjugate; }
};

/** @brief Sphere-to-plane projection of a composed effect's Stage::Project. */
enum class ProjectionKind : uint8_t {
  STEREOGRAPHIC,
  GNOMONIC_FOLDED,
  EQUIRECTANGULAR,
  FOLDED_SINUSOIDAL
};

/** @brief Whether @p projection reads the central-meridian field.
    @param projection Projection kind.
    @return True for the equirectangular and folded-sinusoidal projections. */
constexpr bool uses_central_meridian(ProjectionKind projection) {
  return projection == ProjectionKind::EQUIRECTANGULAR ||
         projection == ProjectionKind::FOLDED_SINUSOIDAL;
}

/** @brief Whether @p projection reads the singularity-fade field. Folded
    sinusoidal has no singular locus and returns fixed weights.
    @param projection Projection kind.
    @return False only for folded sinusoidal. */
constexpr bool uses_singularity_fade(ProjectionKind projection) {
  return projection != ProjectionKind::FOLDED_SINUSOIDAL;
}

/** @brief Optional transfer curve an effect's material stage composes. */
enum class TransferKind : uint8_t { NONE, ISO_CONTOUR };

/** @brief Optional value-dependent coverage stage after sampling. */
enum class FieldCoverageKind : uint8_t { NONE, VALUE_CUTOUT };
