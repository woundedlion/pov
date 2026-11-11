/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/surface.h"
#include "core/render/pullback/authored_lens.h"

/** @file tangent_surface.h
 * @brief Polynomial tangent-field and authored noise displacement policies. */

namespace Pullback::Surface {

/// Catalog noise-basis spellings.
inline constexpr auto &AUTHORED_NOISE_BASIS_IDS = math::NOISE_BASIS_IDS;
/// Fixed integration scheme spellings.
inline constexpr const char *AUTHORED_INTEGRATOR_IDS[] = {"euler", "midpoint",
                                                          "midpoint-2x"};
/// Tangent-length limiter spellings.
inline constexpr const char *AUTHORED_LIMIT_IDS[] = {"smooth", "clamp"};
/// Shared symmetry-group spellings.
inline constexpr auto &AUTHORED_SYMMETRY_IDS =
    Lens::KaleidoscopeParams::SYMMETRY_IDS;
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
  /// Noise basis and integration selectors.
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
  /// Noise basis selector.
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

static_assert(std::size(AUTHORED_INTEGRATOR_IDS) ==
              static_cast<size_t>(Integrator::MIDPOINT_2X) + 1);
static_assert(std::size(AUTHORED_LIMIT_IDS) ==
              static_cast<size_t>(math::TangentLimit::CLAMP) + 1);
static_assert(std::size(AUTHORED_ROSSBY_ORDER_IDS) + 1 <
              std::size(HARMONIC_ORDER_IDS));
static_assert(std::size(POTENTIAL_TERM_IDS) ==
              static_cast<size_t>(math::tangent_field::TermType::DODECAHEDRAL) +
                  1);
static_assert(std::size(POTENTIAL_KIND_IDS) ==
              static_cast<size_t>(math::tangent_field::Kind::POLOIDAL) + 1);
/// Schema fields per scalar-potential term.
inline constexpr size_t POTENTIAL_TERM_FIELDS = 6;
/// Maximum absolute authored term amplitude.
inline constexpr float POTENTIAL_AMPLITUDE_LIMIT = 4.0f;

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
  static constexpr float PHASE_PERIOD = math::TWO_PI_F; ///< Azimuthal period.
  float strength = 0.05f;           ///< Signed fixed-step scale.
  float orientation_azimuth = 0.0f; ///< Polar-axis azimuth.
  float orientation_polar = 0.0f;   ///< Polar-axis colatitude.
  float drift = 0.0f;               ///< Azimuthal advance per frame, radians.
  uint8_t integrator = 0;           ///< Integration scheme.
  std::array<PotentialTermParams, 4> terms{}; ///< Authored terms.
  /// Potential field controls.
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
  /// Integration selector.
  static constexpr auto TOPOLOGY = std::array{TopologyField<PotentialParams>{
      "integrator", &PotentialParams::integrator, AUTHORED_INTEGRATOR_IDS, 0}};
  /// Nested term fields and their mode-only activation gates.
  static constexpr auto EXTRA_SCHEMA = [] {
    std::array<ParamFieldInfo,
               std::tuple_size_v<decltype(terms)> * POTENTIAL_TERM_FIELDS>
        fields{};
    constexpr const char *IDS[][POTENTIAL_TERM_FIELDS] = {
        {"term-1-type", "term-1-kind", "term-1-l", "term-1-m",
         "term-1-amplitude", "term-1-phase"},
        {"term-2-type", "term-2-kind", "term-2-l", "term-2-m",
         "term-2-amplitude", "term-2-phase"},
        {"term-3-type", "term-3-kind", "term-3-l", "term-3-m",
         "term-3-amplitude", "term-3-phase"},
        {"term-4-type", "term-4-kind", "term-4-l", "term-4-m",
         "term-4-amplitude", "term-4-phase"}};
    static_assert(std::size(IDS) == std::tuple_size_v<decltype(terms)>);
    const auto enum_field = []<size_t Count>(const char *id, const char *label,
                                             const char *const(&options)[Count],
                                             uint8_t fallback,
                                             const char *gate = nullptr) {
      static_assert(Count <= UINT8_MAX);
      return ParamFieldInfo{.id = id,
                            .display_name = label,
                            .min = 0.0f,
                            .max = 0.0f,
                            .def = 0.0f,
                            .curve = FieldCurve::SNAP,
                            .topology = true,
                            .enum_count = Count,
                            .enum_def = fallback,
                            .enum_ids = options,
                            .gated_by = gate,
                            .gate_values = gate ? live_values(0) : uint16_t{0}};
    };
    for (size_t i = 0; i < std::size(IDS); ++i) {
      const auto &ids = IDS[i];
      const size_t base = POTENTIAL_TERM_FIELDS * i;
      fields[base] = enum_field(ids[0], "Term Type", POTENTIAL_TERM_IDS, 0);
      fields[base + 1] = enum_field(ids[1], "Kind", POTENTIAL_KIND_IDS, 0);
      fields[base + 2] =
          enum_field(ids[2], "Degree", HARMONIC_ORDER_IDS, 1, ids[0]);
      fields[base + 3] =
          enum_field(ids[3], "Order", HARMONIC_ORDER_IDS, 0, ids[0]);
      fields[base + 4] = {.id = ids[4],
                          .display_name = "Amplitude",
                          .min = -POTENTIAL_AMPLITUDE_LIMIT,
                          .max = POTENTIAL_AMPLITUDE_LIMIT,
                          .def = 0.0f,
                          .curve = FieldCurve::SNAP,
                          .topology = false,
                          .enum_count = 0,
                          .enum_def = 0,
                          .enum_ids = nullptr,
                          .gated_by = nullptr,
                          .gate_values = 0};
      fields[base + 5] = {.id = ids[5],
                          .display_name = "Phase",
                          .min = -math::PI_F,
                          .max = math::PI_F,
                          .def = 0.0f,
                          .curve = FieldCurve::SHORTEST_PERIODIC,
                          .topology = false,
                          .enum_count = 0,
                          .enum_def = 0,
                          .enum_ids = nullptr,
                          .gated_by = ids[0],
                          .gate_values = live_values(0)};
    }
    return fields;
  }();

  /** @brief Addresses one nested term field.
   * @param index Extra-schema field index.
   * @return Parameter storage. */
  HS_FLASH_MEMBER void *extra_address(size_t index) {
    auto &params = *this;
    auto &term = params.terms[index / POTENTIAL_TERM_FIELDS];
    switch (index % POTENTIAL_TERM_FIELDS) {
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
    for (size_t i = 0; i < out.terms.size(); ++i) {
      out.terms[i] = fraction < 1.0f ? from.terms[i] : to.terms[i];
      if (from.terms[i].type == to.terms[i].type &&
          from.terms[i].l == to.terms[i].l &&
          from.terms[i].m == to.terms[i].m) {
        out.terms[i].phase = Control::apply_curve(
            Control::FieldCurve::SHORTEST_PERIODIC, from.terms[i].phase,
            to.terms[i].phase, fraction, -math::PI_F);
      }
    }
  }
};

/** @brief Rossby-Haurwitz field controls. */
struct RossbyHaurwitzParams {
  static constexpr float PHASE_PERIOD = math::TWO_PI_F; ///< Azimuthal period.
  float strength = 0.05f; ///< Signed fixed-step scale.
  float omega = 0.0f;     ///< Zonal peak-speed amplitude.
  float wave = 1.0f;      ///< Wave peak-speed amplitude.
  uint8_t order = 3;      ///< Wave order minus one.
  float drift = 0.0f;     ///< Azimuthal advance per frame, radians.
  uint8_t integrator = 0; ///< Integration scheme.
  /// Rossby–Haurwitz wave controls.
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
  /// Integration selector.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<RossbyHaurwitzParams>{"integrator",
                                          &RossbyHaurwitzParams::integrator,
                                          AUTHORED_INTEGRATOR_IDS, 0},
      TopologyField<RossbyHaurwitzParams>{"order", &RossbyHaurwitzParams::order,
                                          AUTHORED_ROSSBY_ORDER_IDS, 3}};
};

/** @brief Rotating invariant field controls. */
struct SymmetricParams {
  static constexpr float PHASE_PERIOD = math::TWO_PI_F; ///< Rotation period.
  float strength = 0.05f; ///< Signed fixed-step scale.
  float speed = 0.0f;     ///< Field rotation per frame, radians.
  uint8_t integrator = 0; ///< Integration scheme.
  uint8_t symmetry = 0;   ///< Kaleidoscope symmetry index.
  uint8_t kind = 0;       ///< Toroidal or poloidal component.
  /// Symmetric flow controls.
  static constexpr auto FIELDS = std::array{
      Pullback::Field<SymmetricParams>{"strength", &SymmetricParams::strength,
                                       "Strength", -4.0f * math::PI_F,
                                       4.0f * math::PI_F, FieldCurve::SNAP},
      Pullback::Field<SymmetricParams>{"speed", &SymmetricParams::speed,
                                       "Speed", -0.1f, 0.1f, FieldCurve::LERP}};
  /// Integration selector.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<SymmetricParams>{"integrator", &SymmetricParams::integrator,
                                     AUTHORED_INTEGRATOR_IDS, 0},
      TopologyField<SymmetricParams>{"symmetry", &SymmetricParams::symmetry,
                                     AUTHORED_SYMMETRY_IDS, 0},
      TopologyField<SymmetricParams>{"kind", &SymmetricParams::kind,
                                     AUTHORED_FIELD_KIND_IDS, 0}};
};

/** @brief Checks the scalar-potential term domains.
 * @param params Authored field controls.
 * @return Whether all terms are admitted. */
HS_FLASH_MEMBER inline bool valid_terms(const PotentialParams &params) {
  for (const auto &term : params.terms)
    if (term.type >= std::size(POTENTIAL_TERM_IDS) ||
        term.kind >= std::size(POTENTIAL_KIND_IDS) ||
        term.l >= std::size(HARMONIC_ORDER_IDS) ||
        term.m >= std::size(HARMONIC_ORDER_IDS) ||
        (term.type == 0 && term.m > term.l) || !std::isfinite(term.amplitude) ||
        fabsf(term.amplitude) > POTENTIAL_AMPLITUDE_LIMIT ||
        !std::isfinite(term.phase) || fabsf(term.phase) > math::PI_F)
      return false;
  return true;
}

/** @brief Prepares the oriented scalar-potential field.
 * @param params Authored field controls.
 * @param phase Accumulated azimuthal drift, radians.
 * @return Normalized field. */
HS_FLASH_MEMBER inline math::tangent_field::PotentialField
potential_field(const PotentialParams &params, float phase) {
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
  return PotentialField(
      terms, phase,
      math::make_rotation(math::Y_AXIS, math::Vector::from_spherical(
                                            params.orientation_azimuth,
                                            params.orientation_polar)));
}

/** @brief Prepares a Rossby-Haurwitz velocity field.
 * @param omega Zonal peak-speed amplitude.
 * @param wave Wave peak-speed amplitude.
 * @param order Wave order minus one.
 * @param phase Azimuthal drift in radians.
 * @return Normalized field. */
HS_FLASH_MEMBER inline math::tangent_field::PotentialField
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
HS_FLASH_MEMBER inline math::tangent_field::PotentialField
symmetric_field(uint8_t symmetry, uint8_t kind, float phase) {
  using namespace math::tangent_field;
  using Symmetry = Lens::KaleidoscopeSymmetry;
  static_assert(static_cast<uint8_t>(Symmetry::TETRAHEDRAL) ==
                static_cast<uint8_t>(TermType::TETRAHEDRAL));
  static_assert(static_cast<uint8_t>(Symmetry::OCTAHEDRAL) ==
                static_cast<uint8_t>(TermType::OCTAHEDRAL));
  static_assert(static_cast<uint8_t>(Symmetry::DODECAHEDRAL) ==
                static_cast<uint8_t>(TermType::DODECAHEDRAL));
  constexpr auto ORDERS = [] {
    std::array<uint8_t, std::size(Lens::KaleidoscopeParams::SYMMETRY_IDS)>
        orders{};
    orders[static_cast<size_t>(Symmetry::AZIMUTHAL)] = 6;
    orders[static_cast<size_t>(Symmetry::TRIANGULAR_PRISM)] = 3;
    orders[static_cast<size_t>(Symmetry::SQUARE_PRISM)] = 4;
    orders[static_cast<size_t>(Symmetry::PENTAGONAL_PRISM)] = 5;
    orders[static_cast<size_t>(Symmetry::HEXAGONAL_PRISM)] = 6;
    orders[static_cast<size_t>(Symmetry::OCTAGONAL_PRISM)] = 8;
    return orders;
  }();
  static_assert(ORDERS.size() == std::size(AUTHORED_SYMMETRY_IDS));
  std::array<PotentialTerm, 4> terms{};
  const auto type =
      symmetry >= static_cast<uint8_t>(Symmetry::TETRAHEDRAL) &&
              symmetry <= static_cast<uint8_t>(Symmetry::DODECAHEDRAL)
          ? static_cast<TermType>(symmetry)
          : TermType::MODE;
  const uint8_t order = ORDERS[symmetry];
  terms[0] = {type,  static_cast<Kind>(kind),  order,
              order, order % 2 ? -1.0f : 1.0f, 0};
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
  return step_integrated(input, field, strength,
                         static_cast<Integrator>(integrator), track_path);
}

/** @brief Polynomial tangent-field family. */
enum class TangentFieldKind : uint8_t {
  ROSSBY_HAURWITZ, ///< Zonal and travelling-wave field.
  SYMMETRIC,       ///< Symmetry-invariant field.
  POTENTIAL        ///< Authored scalar-potential terms.
};

/** @brief Fixed-step field policy shared by polynomial displacements.
 * @tparam State Parameter and clock provider.
 * @tparam KindV Field kind: Rossby-Haurwitz, symmetric or arbitrary potential. */
template <typename State, TangentFieldKind KindV>
struct TangentField : ApproximationDefaults {
  using Binding = typename State::Binding;              ///< Pipeline binding.
  using FrameState = typename State::FrameState;        ///< Bound frame.
  using Prepared = math::tangent_field::PotentialField; ///< Velocity field.
  /** @brief Prepares this frame's velocity field.
   * @param frame Bound frame.
   * @return Normalized field. */
  static Prepared prepare(const FrameState &frame) {
    const auto &params = State::params(frame);
    if constexpr (KindV == TangentFieldKind::SYMMETRIC)
      return symmetric_field(params.symmetry, params.kind, State::phase(frame));
    else if constexpr (KindV == TangentFieldKind::ROSSBY_HAURWITZ)
      return rossby_haurwitz_field(params.omega, params.wave, params.order,
                                   State::phase(frame));
    else
      return potential_field(params, State::phase(frame));
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
template <typename State>
using RossbyHaurwitz = TangentField<State, TangentFieldKind::ROSSBY_HAURWITZ>;
/** @brief Rotating invariant displacement. @tparam State Field provider. */
template <typename State>
using Symmetric = TangentField<State, TangentFieldKind::SYMMETRIC>;
/** @brief Arbitrary potential displacement. @tparam State Field provider. */
template <typename State>
using Potential = TangentField<State, TangentFieldKind::POTENTIAL>;

} // namespace Pullback::Surface
