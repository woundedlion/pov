/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_parameters.h
 * @brief Parameter-kind enums, uses_* predicates, interpolation, range
 * validation and runtime storage for composed effects.
 */

/**
 * @brief Interpolates one parameter family across a preset transition.
 * @details Driven by the family's field table; each field moves on the curve
 * its descriptor names. Modulation interpolates within a topology and takes
 * the target immediately when topology changes; other untabled members snap
 * to @p b once progress reaches 1.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param t Progress fraction.
 * @return The interpolated family.
 */
template <Pullback::HasFields T>
inline T interpolate(const T &a, const T &b, float t) {
  T out = Pullback::Fields::interpolate(a, b, t);
  if constexpr (requires { out.modulation; })
    out.modulation = Modulation::interpolate(a.modulation, b.modulation, t);
  return out;
}

/** @brief Interpolates signed field angles while retaining potential topology.
 * @param a Outgoing field controls.
 * @param b Incoming field controls.
 * @param t Transition progress.
 * @return Interpolated controls in the signed angle domain. */
inline PotentialParams interpolate(const PotentialParams &a,
                                   const PotentialParams &b, float t) {
  PotentialParams out = Fields::interpolate(a, b, t);
  PotentialParams::sample_extra(out, a, b, t);
  return out;
}

/**
 * @brief Writes interpolated parameter families into existing storage.
 * @tparam Resources Parameter resources.
 * @param out Destination parameter set; may alias either endpoint.
 * @param from Parameters at progress 0.
 * @param to Parameters at progress 1.
 * @param progress Progress fraction, already eased by the caller.
 */
template <typename... Resources>
inline void interpolate_into(
    ComposedDetail::ParameterSet<ComposedDetail::ResourceList<Resources...>>
        &out,
    const ComposedDetail::ParameterSet<
        ComposedDetail::ResourceList<Resources...>> &from,
    const ComposedDetail::ParameterSet<
        ComposedDetail::ResourceList<Resources...>> &to,
    float progress) {
  ((out.template get<Resources::KEY>() =
        interpolate(from.template get<Resources::KEY>(),
                    to.template get<Resources::KEY>(), progress)),
   ...);
}

/** @brief Interpolates a parameter set family by family.
 * @tparam Resources Parameter resources.
 * @param from Outgoing parameters.
 * @param to Incoming parameters.
 * @param progress Eased transition progress.
 * @return Interpolated parameter set. */
template <typename... Resources>
inline ComposedDetail::ParameterSet<ComposedDetail::ResourceList<Resources...>>
interpolate(const ComposedDetail::ParameterSet<
                ComposedDetail::ResourceList<Resources...>> &from,
            const ComposedDetail::ParameterSet<
                ComposedDetail::ResourceList<Resources...>> &to,
            float progress) {
  ComposedDetail::ParameterSet<ComposedDetail::ResourceList<Resources...>> out;
  interpolate_into(out, from, to, progress);
  return out;
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
  if constexpr (requires { value.modulation; })
    return Modulation::admitted(value.modulation) &&
           (T::DYNAMIC || (value.modulation.source == T::SOURCE &&
                           value.modulation.evaluation == T::EVALUATION));
  else
    return Pullback::Fields::valid(value);
}

/** @brief Validates authored numeric and topology controls.
 * @param value Authored controls.
 * @return Whether every control is admitted. */
inline bool valid(const Source::AuthoredProjectedNoiseSourceParams &value) {
  return Pullback::Fields::valid(value) && value.basis <= 2;
}
/** @brief Validates authored numeric and topology controls.
 * @param value Authored controls.
 * @return Whether every control is admitted. */
inline bool valid(const Source::AuthoredSphericalNoiseSourceParams &value) {
  return Pullback::Fields::valid(value) && value.basis <= 2;
}
/** @brief Validates authored numeric and topology controls.
 * @param value Authored controls.
 * @return Whether every control is admitted. */
inline bool valid(const Source::AuthoredTessellationSourceParams &value) {
  return Pullback::Fields::valid(value) && value.kind <= 2;
}
/** @brief Validates authored numeric and topology controls.
 * @param value Authored controls.
 * @return Whether every control is admitted. */
inline bool valid(const Warp::AuthoredWaveShearParams &value) {
  return Pullback::Fields::valid(value) && value.envelope <= 2;
}
/** @brief Validates authored numeric and topology controls.
 * @param value Authored controls.
 * @return Whether every control is admitted. */
inline bool valid(const Warp::AuthoredVectorNoiseParams &value) {
  return Pullback::Fields::valid(value) && value.basis <= 2 &&
         value.envelope <= 2;
}
/** @brief Validates authored numeric and topology controls.
 * @param value Authored controls.
 * @return Whether every control is admitted. */
inline bool valid(const Warp::AuthoredPolarParams &value) {
  return Pullback::Fields::valid(value) && value.mode <= 1 &&
         value.harmonic < Warp::MAX_POLAR_HARMONIC;
}
/** @brief Validates authored numeric and topology controls.
 * @param value Authored controls.
 * @return Whether every control is admitted. */
inline bool valid(const Warp::AuthoredCurlFlowParams &value) {
  return Pullback::Fields::valid(value) && value.basis <= 2 &&
         value.integrator <= 2;
}

/** @brief Validates runtime source and crossing controls.
 * @tparam Family Source-specific controls.
 * @param value Authored source crossing.
 * @return Whether every numeric and topology field is admitted. */
template <typename Family>
inline bool valid(const Source::AuthoredSampleParams<Family> &value) {
  return valid(static_cast<const Family &>(value)) && Fields::valid(value) &&
         value.weight_mode <= 1 &&
         value.coverage_mode <=
             static_cast<uint8_t>(ProjectionCoverageMode::EDGE_FADE_SQUARED);
}

/** @brief Validates ripple controls and scheduler topology.
 * @param value Ripple controls.
 * @return Whether the numeric and topology fields are admitted. */
inline bool valid(const SphereAnimation::BurstRippleParams &value) {
  return Pullback::Fields::valid(value) && value.burst_size <= 3 &&
         value.random_center <= 1;
}

/** @brief Validates Rossby-Haurwitz controls and integration topology.
 * @param value Field controls.
 * @return Whether the controls are admitted. */
inline bool valid(const RossbyHaurwitzParams &value) {
  return Fields::valid(value) && value.order < 7 && value.integrator < 3;
}
/** @brief Validates invariant-field controls and topology.
 * @param value Field controls.
 * @return Whether the controls are admitted. */
inline bool valid(const SymmetricParams &value) {
  return Fields::valid(value) && value.symmetry < 9 && value.kind < 2 &&
         value.integrator < 3;
}
/** @brief Validates arbitrary-potential controls and harmonic domains.
 * @param value Field controls.
 * @return Whether the controls are admitted. */
inline bool valid(const PotentialParams &value) {
  if (!Fields::valid(value) || value.integrator >= 3)
    return false;
  for (const auto &term : value.terms)
    if (term.type >= 4 || term.kind >= 2 || term.l > 8 || term.m > 8 ||
        (term.type == 0 && term.m > term.l) || !std::isfinite(term.amplitude) ||
        fabsf(term.amplitude) > 4 || !std::isfinite(term.phase) ||
        fabsf(term.phase) > math::PI_F)
      return false;
  return true;
}

/** @brief Validates curl-noise numerical topology.
 * @param value Surface controls.
 * @return Whether all fields and topology indices are admitted. */
inline bool valid(const AuthoredSurfaceNoiseParams &value) {
  return Fields::valid(value) && value.basis < 3 && value.integrator < 3 &&
         value.limit < 2;
}
/** @brief Validates directed-noise numerical topology.
 * @param value Surface controls.
 * @return Whether all fields and the noise basis are admitted. */
inline bool valid(const AuthoredDirectSurfaceParams &value) {
  return Fields::valid(value) && value.basis < 3;
}
/** @brief Validates kaleidoscope topology.
 * @param value Lens controls.
 * @return Whether the symmetry is admitted. */
inline bool valid(const Lens::KaleidoscopeParams &value) {
  return value.symmetry < 9;
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

/** @brief Validates authored source controls and selectors.
 * @param p Authored color family.
 * @return Whether all controls and selector ordinals are valid. */
inline bool valid(const Color::AuthoredColorParams &p) {
  return valid(static_cast<const ColorParams &>(p)) && Fields::valid(p) &&
         p.palette_type <=
             static_cast<uint8_t>(AuthoredPaletteType::MUTATING) &&
         p.palette_mode < 3 &&
         p.palette_harmony <= static_cast<uint8_t>(PaletteHarmony::SQUARE) &&
         p.palette_lightness_curve <= static_cast<uint8_t>(AxisCurve::CUP) &&
         p.palette_name < std::size(AUTHORED_PALETTE_NAME_IDS) &&
         p.palette_name_b < std::size(AUTHORED_PALETTE_NAME_IDS) &&
         p.brightness_envelope <= Color::BrightnessEnvelope::DESCENDING;
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

/** @brief Hue-rotation LUT storage; empty when the effect never rotates hue. */
template <bool Enabled> struct OptionalHueRotationLut {};
/** @brief Persistent dodecahedral matrices, empty when disabled.
 * @tparam Enabled Whether reflection matrices are owned. */
template <bool Enabled> struct OptionalDodecahedralMatrices {};
/** @brief Persistent face classifiers, empty when disabled.
 * @tparam Enabled Whether classifier maps are owned. */
template <bool Enabled> struct OptionalDodecahedralFaces {};
/** @brief Immutable cube-face classifiers owned by one effect instance. */
template <> struct OptionalDodecahedralFaces<true> {
  std::array<uint8_t, 6 * lenses::DodecahedralTable::FACE_STEPS *
                          lenses::DodecahedralTable::FACE_STEPS>
      dodecahedral_faces; ///< Face-major exact classifier bytes.
};
/** @brief Immutable reflection matrices owned by one effect instance. */
template <> struct OptionalDodecahedralMatrices<true> {
  std::array<lenses::DodecahedralTable::Matrix,
             lenses::DodecahedralTable::GROUP.count>
      dodecahedral_matrices; ///< Exact copy of the canonical group matrices.
};
/** @brief Interval base storage, empty unless the spec opts in.
 * @tparam Enabled Whether the effect owns palette base caches.
 * @tparam Values Palette-coordinate nodes per cached time row. */
template <bool Enabled, int Values = 64> struct OptionalPaletteLabLut {};
/** @brief Final hue-cache storage, empty unless the spec opts in.
 * @tparam Enabled Whether interval RGB caches are owned.
 * @tparam MidpointOnly Whether only one qualified time row is stored. */
template <bool Enabled, bool MidpointOnly = false>
struct OptionalPaletteRgbLut {};
/** @brief Retained and current final hue tables. */
template <bool MidpointOnly> struct OptionalPaletteRgbLut<true, MidpointOnly> {
  std::array<PaletteRgbHueLutImpl<MidpointOnly>, 2>
      palette_rgb_luts; ///< Ping-pong interval tables.
};
/** @brief Ping-pong interval bases retained across consecutive frames. */
template <int Values> struct OptionalPaletteLabLut<true, Values> {
  std::array<PaletteLabLutStorage<Values>, 2>
      palette_lab_luts;         ///< Previous and current bases.
  uint8_t palette_lab_slot = 0; ///< Most recently baked base.
};
/** @brief Enabled hue-rotation LUT storage. */
template <> struct OptionalHueRotationLut<true> {
  /// Palette colours by palette coordinate and hue-rotation step.
  std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE> hue_rotation_lut;
  /** Palette bake the resident table was built from; 0 matches no bake,
      forcing the first build. */
  uint32_t hue_rotation_lut_bake = 0;
};

/** @brief Projection-walk noise storage; empty when disabled. */
template <bool Enabled> struct ProjectionWalkNoise {};
/** @brief Enabled projection-walk noise storage. */
template <> struct ProjectionWalkNoise<true> {
  FastNoiseLite projection_walk_noise; ///< Noise driving the projection walk.
};

/** @brief Shared outer-camera noise storage; empty when unused. */
template <bool Enabled> struct OuterCameraNoise {};
template <> struct OuterCameraNoise<true> {
  FastNoiseLite outer_walk_noise; ///< Noise driving the camera walk.
};

/** @brief Shared outer-camera clocks; empty when unused. */
template <bool Enabled> struct OuterCameraState {
  math::Quaternion camera_conjugate() const { return {}; }
};
template <> struct OuterCameraState<true> {
  math::Orientation<> outer_walk;       ///< Random-walk orientation.
  math::Quaternion outer_walk_previous; ///< Previous walk orientation.
  math::Quaternion outer_wander;        ///< Accumulated scaled walk deltas.
  math::Quaternion outer_conjugate; ///< Inverse spin-times-wander orientation.
  float camera_spin = 0.0f;         ///< Accumulated Y spin, radians.

  math::Quaternion camera_conjugate() const { return outer_conjugate; }
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
  FOLDED_SINUSOIDAL,
  S2,
  PEIRCE_SQUARE_FAST,
  PEIRCE,
  BONNE,
  AIROCEAN,
  GNOMONIC,
  TETRAHEDRAL,
  CUBIC,
  OCTAHEDRAL,
  DODECAHEDRAL,
  ICOSAHEDRAL,
  NONE
};

/** @brief Whether a projection is a polyhedral atlas.
 * @param projection Projection kind.
 * @return Whether the kind tiles the faces of a Platonic solid. */
constexpr bool is_polyhedral_atlas(ProjectionKind projection) {
  return projection == ProjectionKind::TETRAHEDRAL ||
         projection == ProjectionKind::CUBIC ||
         projection == ProjectionKind::OCTAHEDRAL ||
         projection == ProjectionKind::DODECAHEDRAL ||
         projection == ProjectionKind::ICOSAHEDRAL;
}

/** @brief Whether a projection reads the central meridian.
 * @param projection Projection kind.
 * @return Whether meridian controls affect the selected kernel. */
constexpr bool uses_central_meridian(ProjectionKind projection) {
  return projection == ProjectionKind::EQUIRECTANGULAR ||
         projection == ProjectionKind::FOLDED_SINUSOIDAL ||
         projection == ProjectionKind::S2 ||
         projection == ProjectionKind::PEIRCE ||
         projection == ProjectionKind::BONNE ||
         projection == ProjectionKind::AIROCEAN ||
         is_polyhedral_atlas(projection);
}

/** @brief Whether a projection reads singularity attenuation.
 * @param projection Projection kind.
 * @return Whether singularity controls affect the selected kernel. */
constexpr bool uses_singularity_fade(ProjectionKind projection) {
  return projection == ProjectionKind::STEREOGRAPHIC ||
         projection == ProjectionKind::GNOMONIC_FOLDED ||
         projection == ProjectionKind::GNOMONIC ||
         projection == ProjectionKind::PEIRCE ||
         projection == ProjectionKind::EQUIRECTANGULAR ||
         projection == ProjectionKind::PEIRCE_SQUARE_FAST;
}

/** @brief Optional transfer curve an effect's material stage composes. */
enum class TransferKind : uint8_t { NONE, ISO_CONTOUR, RIDGE, SMOOTH_BANDS };

/** @brief Optional value-dependent coverage stage after sampling. */
enum class FieldCoverageKind : uint8_t { NONE, VALUE_CUTOUT };

/** @brief Extended controls for atlas and selectable-hemisphere projections. */
struct AtlasProjectionParams : ProjectionParams {
  float coordinate_scale = 1.0f;               ///< Plane coordinate scale.
  float layout_scroll = 0.0f;                  ///< Peirce strip translation.
  float standard_parallel = math::PI_F / 4.0f; ///< Bonne reference latitude.
  float area_correction = 0.0f; ///< Polyhedral atlas blend toward equal area.
  uint8_t layout = 1;           ///< Peirce or Airocean layout.
  uint8_t hemisphere = 0;       ///< Gnomonic or Bonne hemisphere.
  uint8_t frame = 1;            ///< Identity or animated projection frame.
  static constexpr auto FIELDS = concat_fields<AtlasProjectionParams>(
      ProjectionParams::FIELDS,
      std::array{
          Field<AtlasProjectionParams>{
              "coordinate-scale", &AtlasProjectionParams::coordinate_scale,
              "Coordinate Scale", 0.25f, 4.0f, FieldCurve::LOG_POSITIVE,
              FieldGate::COORDINATE_SCALE},
          Field<AtlasProjectionParams>{
              "layout-scroll", &AtlasProjectionParams::layout_scroll,
              "Layout Scroll", -1.0f, 1.0f, FieldCurve::LERP,
              FieldGate::LAYOUT_SCROLL},
          Field<AtlasProjectionParams>{
              "standard-parallel", &AtlasProjectionParams::standard_parallel,
              "Standard Parallel", 0.0f, math::PI_F / 2.0f, FieldCurve::LERP,
              FieldGate::STANDARD_PARALLEL},
          Field<AtlasProjectionParams>{
              "area-correction", &AtlasProjectionParams::area_correction,
              "Area Correction", 0.0f, 1.0f, FieldCurve::LERP,
              FieldGate::AREA_CORRECTION}});
  static constexpr const char *FRAME_IDS[] = {"identity", "spin-wander"};
  static constexpr const char *LAYOUT_IDS[] = {"diamond", "square",
                                               "horizontal", "vertical"};
  static constexpr const char *HEMISPHERE_IDS[] = {"folded", "front", "back"};
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AtlasProjectionParams>{
          "frame", &AtlasProjectionParams::frame, FRAME_IDS, 1},
      TopologyField<AtlasProjectionParams>{
          "layout", &AtlasProjectionParams::layout, LAYOUT_IDS, 1},
      TopologyField<AtlasProjectionParams>{
          "hemisphere", &AtlasProjectionParams::hemisphere, HEMISPHERE_IDS, 0}};
};
