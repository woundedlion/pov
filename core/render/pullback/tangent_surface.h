/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/surface.h"

/** @file tangent_surface.h
 * @brief Polynomial tangent-field displacement policies. */

namespace Pullback::Surface {

/// Catalog noise-basis spellings.
inline constexpr const char *AUTHORED_NOISE_BASIS_IDS[] = {"simplex", "fbm3",
                                                           "ridged3"};
/// Fixed integration scheme spellings.
inline constexpr const char *AUTHORED_INTEGRATOR_IDS[] = {"euler", "midpoint",
                                                          "midpoint-2x"};
/// Tangent-length limiter spellings.
inline constexpr const char *AUTHORED_LIMIT_IDS[] = {"smooth", "clamp"};
/// Shared symmetry-group spellings.
inline constexpr const char *AUTHORED_SYMMETRY_IDS[] = {
    "azimuthal",        "tetrahedral",      "octahedral",
    "dodecahedral",     "triangular-prism", "square-prism",
    "pentagonal-prism", "hexagonal-prism",  "octagonal-prism"};
/// Helmholtz component spellings.
inline constexpr const char *AUTHORED_FIELD_KIND_IDS[] = {"toroidal",
                                                          "poloidal"};
/// Rossby wave-order spellings.
inline constexpr const char *AUTHORED_ROSSBY_ORDER_IDS[] = {
    "one", "two", "three", "four", "five", "six", "seven"};

/** @brief Curl-noise controls with authored numerical topology. */
struct AuthoredSurfaceNoiseParams : SurfaceNoiseParams {
  uint8_t basis =
      static_cast<uint8_t>(math::NoiseBasis::SIMPLEX); ///< Noise basis.
  uint8_t integrator = 0;                              ///< Integration scheme.
  uint8_t limit = 0; ///< Tangent-length limit.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AuthoredSurfaceNoiseParams>{
          "basis", &AuthoredSurfaceNoiseParams::basis, AUTHORED_NOISE_BASIS_IDS,
          0},
      TopologyField<AuthoredSurfaceNoiseParams>{
          "integrator", &AuthoredSurfaceNoiseParams::integrator,
          AUTHORED_INTEGRATOR_IDS, 0},
      TopologyField<AuthoredSurfaceNoiseParams>{
          "limit", &AuthoredSurfaceNoiseParams::limit, AUTHORED_LIMIT_IDS, 0}};
};

/** @brief Directed-noise controls with an authored noise basis. */
struct AuthoredDirectSurfaceParams : DirectSurfaceParams {
  uint8_t basis =
      static_cast<uint8_t>(math::NoiseBasis::SIMPLEX); ///< Noise basis.
  static constexpr auto TOPOLOGY =
      std::array{TopologyField<AuthoredDirectSurfaceParams>{
          "basis", &AuthoredDirectSurfaceParams::basis,
          AUTHORED_NOISE_BASIS_IDS, 0}};
};

/** @brief Curl-noise policy with runtime numerical topology.
 * @tparam State Surface provider. */
template <typename State> struct ConfiguredCurl : PrepareFromState<State> {
  using Binding = typename State::Binding;       ///< Pipeline binding.
  using FrameState = typename State::FrameState; ///< Bound frame.
  using Prepared = StatePrepared<State>;         ///< Noise-loop preparation.
  /** @brief Displaces a direction through the configured curl field.
   * @param input Unit direction.
   * @param frame Bound frame.
   * @param prepared Prepared noise loop.
   * @return Displaced direction and angular travel. */
  static SurfaceResult apply(const math::Vector &input, const FrameState &frame,
                             const Prepared &prepared) {
    const auto &params = State::params(frame);
    return curl_noise(input, State::noise(frame),
                      static_cast<math::NoiseBasis>(params.basis),
                      static_cast<Integrator>(params.integrator),
                      static_cast<math::TangentLimit>(params.limit),
                      params.scale, prepared.loop_offset, params.strength,
                      State::path_length_required(frame));
  }
};

/** @brief Directed-noise policy with runtime basis.
 * @tparam State Surface provider. */
template <typename State> struct ConfiguredDirect : PrepareFromState<State> {
  using Binding = typename State::Binding;       ///< Pipeline binding.
  using FrameState = typename State::FrameState; ///< Bound frame.
  using Prepared = StatePrepared<State>;         ///< Noise-loop preparation.
  /** @brief Displaces a direction through the configured noise field.
   * @param input Unit direction.
   * @param frame Bound frame.
   * @param prepared Prepared noise loop and steering frame.
   * @return Displaced direction and angular travel. */
  static SurfaceResult apply(const math::Vector &input, const FrameState &frame,
                             const Prepared &prepared) {
    const auto &params = State::params(frame);
    return direct_noise(input, State::noise(frame),
                        static_cast<math::NoiseBasis>(params.basis),
                        params.scale, prepared.loop_offset, params.strength,
                        prepared.direction_cos, prepared.direction_sin,
                        State::path_length_required(frame));
  }
};

/// Scalar-potential topology spellings.
inline constexpr const char *POTENTIAL_TERM_IDS[] = {
    "mode", "tetrahedral", "octahedral", "dodecahedral"};
/// Potential Helmholtz component spellings.
inline constexpr auto &POTENTIAL_KIND_IDS = AUTHORED_FIELD_KIND_IDS;
/// Harmonic degree and order spellings.
inline constexpr const char *HARMONIC_ORDER_IDS[] = {
    "zero", "one", "two", "three", "four", "five", "six", "seven", "eight"};

/** @brief One authored scalar-potential term. */
struct PotentialTermParams {
  uint8_t type = 0;       ///< Scalar potential topology.
  uint8_t kind = 0;       ///< Toroidal or poloidal component.
  uint8_t l = 1;          ///< Harmonic degree.
  uint8_t m = 0;          ///< Harmonic order.
  float amplitude = 0.0f; ///< Peak-speed amplitude.
  float phase = 0.0f;     ///< Mode phase, radians.
};

/** @brief Four-term polynomial field controls. */
struct PotentialParams {
  float strength = 0.05f;           ///< Signed fixed-step scale.
  float orientation_azimuth = 0.0f; ///< Polar-axis azimuth.
  float orientation_polar = 0.0f;   ///< Polar-axis colatitude.
  float drift = 0.0f;               ///< Azimuthal advance per frame, radians.
  uint8_t integrator = 0;           ///< Integration scheme.
  std::array<PotentialTermParams, 4> terms{}; ///< Authored terms.
  static constexpr auto FIELDS = std::array{
      Pullback::Field<PotentialParams>{"strength", &PotentialParams::strength,
                                       "Strength", -math::PI_F, math::PI_F,
                                       FieldCurve::SNAP},
      Pullback::Field<PotentialParams>{
          "orientation-azimuth", &PotentialParams::orientation_azimuth,
          "Orientation Azimuth", -math::PI_F, math::PI_F,
          FieldCurve::SHORTEST_PERIODIC},
      Pullback::Field<PotentialParams>{
          "orientation-polar", &PotentialParams::orientation_polar,
          "Orientation Polar", 0.0f, math::PI_F, FieldCurve::LERP},
      Pullback::Field<PotentialParams>{"drift", &PotentialParams::drift,
                                       "Drift", -0.1f, 0.1f, FieldCurve::LERP}};
  static constexpr auto TOPOLOGY = std::array{TopologyField<PotentialParams>{
      "integrator", &PotentialParams::integrator, AUTHORED_INTEGRATOR_IDS, 0}};
  /// Nested term fields and their mode-only activation gates.
  static constexpr auto EXTRA_SCHEMA = [] {
    std::array<ParamFieldInfo, 24> fields{};
    constexpr const char *IDS[4][6] = {
        {"term-1-type", "term-1-kind", "term-1-l", "term-1-m",
         "term-1-amplitude", "term-1-phase"},
        {"term-2-type", "term-2-kind", "term-2-l", "term-2-m",
         "term-2-amplitude", "term-2-phase"},
        {"term-3-type", "term-3-kind", "term-3-l", "term-3-m",
         "term-3-amplitude", "term-3-phase"},
        {"term-4-type", "term-4-kind", "term-4-l", "term-4-m",
         "term-4-amplitude", "term-4-phase"}};
    for (size_t i = 0; i < 4; ++i) {
      fields[6 * i] = {IDS[i][0], "Term Type",        0,       0,
                       0,         FieldCurve::SNAP,   true,    4,
                       0,         POTENTIAL_TERM_IDS, nullptr, 0};
      fields[6 * i + 1] = {IDS[i][1],        "Kind", 0, 0, 0,
                           FieldCurve::SNAP, true,   2, 0, POTENTIAL_KIND_IDS,
                           nullptr,          0};
      fields[6 * i + 2] = {
          IDS[i][2], "Degree",           0,         0,
          0,         FieldCurve::SNAP,   true,      9,
          1,         HARMONIC_ORDER_IDS, IDS[i][0], live_values(0)};
      fields[6 * i + 3] = {
          IDS[i][3],        "Order",       0, 0, 0,
          FieldCurve::SNAP, true,          9, 0, HARMONIC_ORDER_IDS,
          IDS[i][0],        live_values(0)};
      fields[6 * i + 4] = {IDS[i][4], "Amplitude",      -4,      4,
                           0,         FieldCurve::SNAP, false,   0,
                           0,         nullptr,          nullptr, 0};
      fields[6 * i + 5] = {IDS[i][5],  "Phase",   -math::PI_F,
                           math::PI_F, 0,         FieldCurve::SHORTEST_PERIODIC,
                           false,      0,         0,
                           nullptr,    IDS[i][0], live_values(0)};
    }
    return fields;
  }();

  /** @brief Addresses one nested term field.
   * @param index Extra-schema field index.
   * @return Parameter storage. */
  HS_COLD_MEMBER void *extra_address(size_t index) {
    auto &params = *this;
    auto &term = params.terms[index / 6];
    switch (index % 6) {
    case 0:
      return &term.type;
    case 1:
      return &term.kind;
    case 2:
      return &term.l;
    case 3:
      return &term.m;
    case 4:
      return &term.amplitude;
    default:
      return &term.phase;
    }
  }

  /** @brief Samples nested mode phases while retaining endpoint topology.
   * @param out Destination controls.
   * @param from Start controls.
   * @param to End controls.
   * @param fraction Frame fraction. */
  static void sample_extra(PotentialParams &out, const PotentialParams &from,
                           const PotentialParams &to, float fraction) {
    if (out.orientation_azimuth > math::PI_F)
      out.orientation_azimuth -= math::TWO_PI_F;
    for (size_t i = 0; i < 4; ++i) {
      out.terms[i] = fraction < 1.0f ? from.terms[i] : to.terms[i];
      if (from.terms[i].type == to.terms[i].type &&
          from.terms[i].l == to.terms[i].l &&
          from.terms[i].m == to.terms[i].m) {
        out.terms[i].phase = Control::apply_curve(
            Control::FieldCurve::SHORTEST_PERIODIC, from.terms[i].phase,
            to.terms[i].phase, fraction);
        if (out.terms[i].phase > math::PI_F)
          out.terms[i].phase -= math::TWO_PI_F;
      }
    }
  }
};

/** @brief Rossby-Haurwitz field controls. */
struct RossbyHaurwitzParams {
  float strength = 0.05f; ///< Signed fixed-step scale.
  float omega = 0.0f;     ///< Zonal peak-speed amplitude.
  float wave = 1.0f;      ///< Wave peak-speed amplitude.
  uint8_t order = 3;      ///< Wave order minus one.
  float drift = 0.0f;     ///< Azimuthal advance per frame, radians.
  uint8_t integrator = 0; ///< Integration scheme.
  static constexpr auto FIELDS =
      std::array{Pullback::Field<RossbyHaurwitzParams>{
                     "strength", &RossbyHaurwitzParams::strength, "Strength",
                     -math::PI_F, math::PI_F, FieldCurve::SNAP},
                 Pullback::Field<RossbyHaurwitzParams>{
                     "omega", &RossbyHaurwitzParams::omega, "Zonal Amplitude",
                     -4, 4, FieldCurve::SNAP},
                 Pullback::Field<RossbyHaurwitzParams>{
                     "wave", &RossbyHaurwitzParams::wave, "Wave Amplitude", -4,
                     4, FieldCurve::SNAP},
                 Pullback::Field<RossbyHaurwitzParams>{
                     "drift", &RossbyHaurwitzParams::drift, "Drift", -0.1f,
                     0.1f, FieldCurve::LERP}};
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<RossbyHaurwitzParams>{"integrator",
                                          &RossbyHaurwitzParams::integrator,
                                          AUTHORED_INTEGRATOR_IDS, 0},
      TopologyField<RossbyHaurwitzParams>{"order", &RossbyHaurwitzParams::order,
                                          AUTHORED_ROSSBY_ORDER_IDS, 3}};
};

/** @brief Rotating invariant field controls. */
struct SymmetricParams {
  float strength = 0.05f; ///< Signed fixed-step scale.
  float speed = 0.0f;     ///< Field rotation per frame, radians.
  uint8_t integrator = 0; ///< Integration scheme.
  uint8_t symmetry = 0;   ///< Kaleidoscope symmetry index.
  uint8_t kind = 0;       ///< Toroidal or poloidal component.
  static constexpr auto FIELDS = std::array{
      Pullback::Field<SymmetricParams>{"strength", &SymmetricParams::strength,
                                       "Strength", -4.0f * math::PI_F,
                                       4.0f * math::PI_F, FieldCurve::SNAP},
      Pullback::Field<SymmetricParams>{"speed", &SymmetricParams::speed,
                                       "Speed", -0.1f, 0.1f, FieldCurve::LERP}};
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<SymmetricParams>{"integrator", &SymmetricParams::integrator,
                                     AUTHORED_INTEGRATOR_IDS, 0},
      TopologyField<SymmetricParams>{"symmetry", &SymmetricParams::symmetry,
                                     AUTHORED_SYMMETRY_IDS, 0},
      TopologyField<SymmetricParams>{"kind", &SymmetricParams::kind,
                                     AUTHORED_FIELD_KIND_IDS, 0}};
};

/** @brief Prepares a Rossby-Haurwitz velocity field.
 * @param omega Zonal peak-speed amplitude.
 * @param wave Wave peak-speed amplitude.
 * @param order Wave order minus one.
 * @param phase Azimuthal drift in radians.
 * @return Normalized field. */
inline math::tangent_field::PotentialField
rossby_haurwitz_field(float omega, float wave, uint8_t order, float phase) {
  using namespace math::tangent_field;
  std::array<PotentialTerm, 4> terms{};
  terms[0] = {TermType::MODE, Kind::TOROIDAL, 1, 0, -omega, 0};
  terms[1] = {TermType::MODE,
              Kind::TOROIDAL,
              static_cast<uint8_t>(order + 2),
              static_cast<uint8_t>(order + 1),
              wave * (order % 2 ? 1.0f : -1.0f),
              0};
  return PotentialField(terms, phase);
}

/** @brief Prepares a rotating invariant velocity field.
 * @param symmetry Kaleidoscope symmetry index.
 * @param kind Toroidal or poloidal component.
 * @param phase Field rotation in radians.
 * @return Normalized field. */
HS_COLD_MEMBER inline math::tangent_field::PotentialField
symmetric_field(uint8_t symmetry, uint8_t kind, float phase) {
  using namespace math::tangent_field;
  constexpr uint8_t ORDERS[] = {6, 0, 0, 0, 3, 4, 5, 6, 8};
  std::array<PotentialTerm, 4> terms{};
  const auto TYPE = symmetry >= 1 && symmetry <= 3
                        ? static_cast<TermType>(symmetry)
                        : TermType::MODE;
  const uint8_t ORDER = ORDERS[symmetry];
  terms[0] = {TYPE,  static_cast<Kind>(kind),  ORDER,
              ORDER, ORDER % 2 ? -1.0f : 1.0f, 0};
  return PotentialField(terms, 0.0f, math::make_rotation(math::Y_AXIS, phase));
}

/** @brief Dispatches an authored fixed-step integration scheme.
 * @param input Unit direction.
 * @param field Prepared velocity field.
 * @param strength Signed displacement scale.
 * @param integrator Integration scheme index.
 * @param track_path Whether to accumulate angular travel.
 * @return Displaced direction and angular travel. */
HS_HOT_FLASH_MEMBER inline SurfaceResult
step_tangent_field(const math::Vector &input,
                   const math::tangent_field::PotentialField &field,
                   float strength, uint8_t integrator, bool track_path) {
  using Field = math::tangent_field::PotentialField;
  switch (static_cast<Integrator>(integrator)) {
  case Integrator::EULER:
    return step_field<Field, Euler>(input, field, strength, track_path);
  case Integrator::MIDPOINT:
    return step_field<Field, Midpoint>(input, field, strength, track_path);
  default:
    return step_field<Field, Midpoint2>(input, field, strength, track_path);
  }
}

/** @brief Fixed-step field policy shared by polynomial displacements.
 * @tparam State Parameter and clock provider.
 * @tparam KindV Field kind: Rossby-Haurwitz, symmetric or arbitrary potential. */
template <typename State, uint8_t KindV>
struct TangentField : ApproximationDefaults {
  using Binding = typename State::Binding;              ///< Pipeline binding.
  using FrameState = typename State::FrameState;        ///< Bound frame.
  using Prepared = math::tangent_field::PotentialField; ///< Velocity field.
  /** @brief Prepares this frame's velocity field.
   * @param frame Bound frame.
   * @return Normalized field. */
  static Prepared prepare(const FrameState &frame) {
    const auto &params = State::params(frame);
    if constexpr (KindV == 1)
      return symmetric_field(params.symmetry, params.kind, State::phase(frame));
    else if constexpr (KindV == 0)
      return rossby_haurwitz_field(params.omega, params.wave, params.order,
                                   State::phase(frame));
    else {
      using namespace math::tangent_field;
      std::array<PotentialTerm, 4> terms{};
      for (size_t i = 0; i < terms.size(); ++i) {
        const auto &term = params.terms[i];
        terms[i] = {static_cast<TermType>(term.type),
                    static_cast<Kind>(term.kind),
                    term.l,
                    term.m,
                    term.amplitude,
                    term.phase};
      }
      return Prepared(
          terms, State::phase(frame),
          math::make_rotation(math::Y_AXIS, math::Vector::from_spherical(
                                                params.orientation_azimuth,
                                                params.orientation_polar)));
    }
  }
  /** @brief Applies the authored integrator to a direction.
   * @param input Unit direction.
   * @param frame Bound frame.
   * @param prepared Prepared velocity field.
   * @return Displaced direction and angular travel. */
  static SurfaceResult apply(const math::Vector &input, const FrameState &frame,
                             const Prepared &prepared) {
    const auto &params = State::params(frame);
    return step_tangent_field(input, prepared, params.strength,
                              params.integrator,
                              State::path_length_required(frame));
  }
};

/** @brief Rossby-Haurwitz displacement. @tparam State Field provider. */
template <typename State> using RossbyHaurwitz = TangentField<State, 0>;
/** @brief Rotating invariant displacement. @tparam State Field provider. */
template <typename State> using Symmetric = TangentField<State, 1>;
/** @brief Arbitrary potential displacement. @tparam State Field provider. */
template <typename State> using Potential = TangentField<State, 2>;

} // namespace Pullback::Surface
