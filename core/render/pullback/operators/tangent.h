/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER
#include "core/math/tangent_field.h"
#include "core/render/pullback/operators/sphere.h"

/** @file tangent.h
 * @brief Fixed-step polynomial tangent-field operators. */

namespace Pullback::Interp::Op {

/** @brief Scalar-potential topology names. */
inline constexpr const char *POTENTIAL_TERM_IDS[] = {
    "mode", "tetrahedral", "octahedral", "dodecahedral"};
/** @brief Helmholtz component topology names. */
inline constexpr const char *POTENTIAL_KIND_IDS[] = {"toroidal", "poloidal"};
/** @brief Harmonic degree and order values, indexed by the integer value. */
inline constexpr const char *HARMONIC_ORDER_IDS[] = {
    "zero", "one", "two", "three", "four", "five", "six", "seven", "eight"};
/** @brief Rossby wave-order values; the stored index represents R = index + 1. */
inline constexpr const char *ROSSBY_ORDER_IDS[] = {
    "one", "two", "three", "four", "five", "six", "seven"};
/** @brief Rejected field configuration diagnostic. */
inline constexpr const char *FIELD_RANGE_WARNING =
    "Field parameters must be finite and strength times peak-speed bound must not exceed pi.";

/** @brief Conservative strength before a fixed-step map can fold.
 * @param field Prepared tangent field.
 * @return Advisory strength, or infinity for a zero field. */
inline float
field_fold_strength(const math::tangent_field::PotentialField &field) {
  return math::tangent_field::fold_strength(field.speed_bound(),
                                            field.lipschitz_bound());
}

/** @brief Accumulated azimuthal drift clock. */
struct FieldPhaseState {
  float phase = 0.0f; ///< Azimuthal phase, radians modulo 2pi.
};

/** @brief Retains field drift endpoints.
 * @param from Previous clock.
 * @param to Current clock.
 * @return Prepared circular interval. */
inline Spacetime::ClockInterval temporal_interval(const FieldPhaseState &from,
                                                  const FieldPhaseState &to) {
  return {from.phase, to.phase, math::TWO_PI_F};
}

/** @brief Samples field drift at column time.
 * @tparam Params Operator controls.
 * @param out Column clock.
 * @param interval Retained endpoints.
 * @param fraction Column fraction.
 */
template <typename Params>
inline void temporal_sample(FieldPhaseState &out,
                            const Spacetime::ClockInterval &interval,
                            const Params &, float fraction) {
  out.phase = interval.sample(fraction);
}

/** @brief Authored scalar-potential term. */
struct PotentialTermParams {
  uint8_t type = 0;       ///< Scalar potential topology.
  uint8_t kind = 0;       ///< Toroidal or poloidal component.
  uint8_t l = 1;          ///< Harmonic degree.
  uint8_t m = 0;          ///< Harmonic order.
  float amplitude = 0.0f; ///< Peak-speed amplitude.
  float phase = 0.0f;     ///< Mode phase, radians.
};

/** @brief Four-term polynomial field controls. */
struct PotentialDisplaceParams {
  float strength = 0.05f;           ///< Signed fixed-step scale.
  float orientation_azimuth = 0.0f; ///< Field polar-axis azimuth.
  float orientation_polar = 0.0f;   ///< Field polar-axis colatitude.
  float drift = 0.0f;               ///< Drift phase advance per frame, radians.
  uint8_t integrator = 0;           ///< Integration scheme.
  std::array<PotentialTermParams, 4> terms{}; ///< Authored term slots.

  static constexpr auto FIELDS = std::array{
      Field<PotentialDisplaceParams>{
          "strength", &PotentialDisplaceParams::strength, "Strength",
          -math::PI_F, math::PI_F, FieldCurve::SNAP},
      Field<PotentialDisplaceParams>{
          "orientation-azimuth", &PotentialDisplaceParams::orientation_azimuth,
          "Orientation Azimuth", -math::PI_F, math::PI_F,
          FieldCurve::SHORTEST_PERIODIC},
      Field<PotentialDisplaceParams>{
          "orientation-polar", &PotentialDisplaceParams::orientation_polar,
          "Orientation Polar", 0.0f, math::PI_F, FieldCurve::LERP},
      Field<PotentialDisplaceParams>{"drift", &PotentialDisplaceParams::drift,
                                     "Drift", -0.1f, 0.1f, FieldCurve::LERP}};
  static constexpr auto TOPOLOGY =
      std::array{TopologyField<PotentialDisplaceParams>{
          "integrator", &PotentialDisplaceParams::integrator,
          SURFACE_INTEGRATOR_IDS, 0}};
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
  void *extra_address(size_t index) {
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
  static void sample_extra(PotentialDisplaceParams &out,
                           const PotentialDisplaceParams &from,
                           const PotentialDisplaceParams &to, float fraction) {
    for (size_t i = 0; i < 4; ++i) {
      out.terms[i] = fraction < 1.0f ? from.terms[i] : to.terms[i];
      if (from.terms[i].type == to.terms[i].type &&
          from.terms[i].l == to.terms[i].l && from.terms[i].m == to.terms[i].m)
        out.terms[i].phase = Control::apply_curve(
            Control::FieldCurve::SHORTEST_PERIODIC, from.terms[i].phase,
            to.terms[i].phase, fraction);
    }
  }
};

/** @brief Rossby-Haurwitz field controls. */
struct RossbyHaurwitzParams {
  float strength = 0.05f; ///< Signed fixed-step scale.
  float omega = 0.0f;     ///< Zonal peak-speed amplitude.
  float wave = 1.0f;      ///< Wave peak-speed amplitude.
  uint8_t order = 3;      ///< Wave order R = order + 1.
  float drift = 0.0f;     ///< Drift phase advance per frame.
  uint8_t integrator = 0; ///< Integration scheme.
  static constexpr auto FIELDS = std::array{
      Field<RossbyHaurwitzParams>{"strength", &RossbyHaurwitzParams::strength,
                                  "Strength", -math::PI_F, math::PI_F,
                                  FieldCurve::SNAP},
      Field<RossbyHaurwitzParams>{"omega", &RossbyHaurwitzParams::omega,
                                  "Zonal Amplitude", -4, 4, FieldCurve::SNAP},
      Field<RossbyHaurwitzParams>{"wave", &RossbyHaurwitzParams::wave,
                                  "Wave Amplitude", -4, 4, FieldCurve::SNAP},
      Field<RossbyHaurwitzParams>{"drift", &RossbyHaurwitzParams::drift,
                                  "Drift", -0.1f, 0.1f, FieldCurve::LERP}};
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<RossbyHaurwitzParams>{"integrator",
                                          &RossbyHaurwitzParams::integrator,
                                          SURFACE_INTEGRATOR_IDS, 0},
      TopologyField<RossbyHaurwitzParams>{"order", &RossbyHaurwitzParams::order,
                                          ROSSBY_ORDER_IDS, 3}};
};

/** @brief Fold-aligned invariant field controls. */
struct SymmetricDisplaceParams {
  float strength = 0.05f; ///< Signed fixed-step scale.
  float amplitude = 1.0f; ///< Peak-speed amplitude.
  uint8_t integrator = 0; ///< Integration scheme.
  uint8_t symmetry = 0;   ///< Shared kaleidoscope symmetry.
  uint8_t kind = 0;       ///< Helmholtz component.
  static constexpr auto FIELDS =
      std::array{Field<SymmetricDisplaceParams>{
                     "strength", &SymmetricDisplaceParams::strength, "Strength",
                     -math::PI_F, math::PI_F, FieldCurve::SNAP},
                 Field<SymmetricDisplaceParams>{
                     "amplitude", &SymmetricDisplaceParams::amplitude,
                     "Amplitude", -4, 4, FieldCurve::SNAP}};
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<SymmetricDisplaceParams>{
          "integrator", &SymmetricDisplaceParams::integrator,
          SURFACE_INTEGRATOR_IDS, 0},
      TopologyField<SymmetricDisplaceParams>{"symmetry",
                                             &SymmetricDisplaceParams::symmetry,
                                             KALEIDOSCOPE_SYMMETRY_IDS, 0},
      TopologyField<SymmetricDisplaceParams>{
          "kind", &SymmetricDisplaceParams::kind, POTENTIAL_KIND_IDS, 0}};
};

/** @brief Shared operator lifecycle over a fixed-step field.
 * @tparam Derived Concrete field operator.
 * @tparam ParamsT Authored controls.
 * @tparam StateT Instance clock state. */
template <typename Derived, typename ParamsT, typename StateT = FieldPhaseState>
struct FieldDisplaceModel : ValueStateModel<StateT> {
  using State = StateT;
  using Params = ParamsT;
  using Input = SphereSample;
  using Output = SphereSample;
  using Prepared = math::tangent_field::PotentialField;
  /// Interpolation admission contract for the fixed-step arc cap.
  static constexpr const char *ADMISSIBILITY_CONVEXITY =
      "Strength and amplitudes snap together; interpolated orientations and phases preserve speed bounds.";
  static void advance(State &state, const Params &params) {
    if constexpr (requires { params.drift; })
      state.phase = fmodf(state.phase + params.drift, math::TWO_PI_F);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    if constexpr (requires { state.phase; })
      return Derived::field(params, state.phase);
    else
      return Derived::field(params, 0.0f);
  }
  /** @brief Reads the conservative fold threshold.
   * @param params Authored field controls.
   * @return Advisory absolute strength, or infinity for a zero field. */
  static float fold_strength(const Params &params) {
    return field_fold_strength(Derived::field(params, 0.0f));
  }
  /** @brief Checks field domains and the uniform arc cap.
   * @param params Authored field controls.
   * @return Null when admitted, otherwise the range diagnostic. */
  static const char *validate(const Params &params) {
    for (const auto &f : Params::FIELDS)
      if (!std::isfinite(params.*f.member) || params.*f.member < f.min ||
          params.*f.member > f.max)
        return FIELD_RANGE_WARNING;
    if (params.integrator > 2 || !Derived::valid_terms(params))
      return FIELD_RANGE_WARNING;
    return fabsf(params.strength) *
                       Derived::field(params, 0.0f).speed_bound() <=
                   math::PI_F
               ? nullptr
               : FIELD_RANGE_WARNING;
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &params, const Prepared &prepared) {
    return Kernel::displace(
        input, step(input.dir, prepared, params.strength, params.integrator));
  }
  /** @brief Dispatches one fixed field step.
   * @param p Unit start direction.
   * @param field Prepared velocity and bounds.
   * @param strength Signed displacement scale.
   * @param integrator Integration scheme index.
   * @return Displaced direction and accumulated arc length. */
  static SurfaceResult step(const math::Vector &p, const Prepared &field,
                            float strength, uint8_t integrator) {
    switch (static_cast<Surface::Integrator>(integrator)) {
    case Surface::Integrator::EULER:
      return Surface::step_field<Prepared, Surface::Euler>(p, field, strength,
                                                           true);
    case Surface::Integrator::MIDPOINT:
      return Surface::step_field<Prepared, Surface::Midpoint>(p, field,
                                                              strength, true);
    default:
      return Surface::step_field<Prepared, Surface::Midpoint2>(p, field,
                                                               strength, true);
    }
  }
};

/** @brief Four scalar-potential terms assembled into one displacement. */
struct DisplacePotential
    : FieldDisplaceModel<DisplacePotential, PotentialDisplaceParams> {
  static constexpr const char *ID = "sphere.displace.potential.v1";
  static constexpr const char *NAME = "Potential Displace";
  /** @brief Checks the term topology and parameter domains.
   * @param params Authored field controls.
   * @return Whether the term configuration is admitted. */
  static bool valid_terms(const Params &params) {
    for (const auto &t : params.terms)
      if (t.type > 3 || t.kind > 1 || !std::isfinite(t.amplitude) ||
          fabsf(t.amplitude) > 4 || !std::isfinite(t.phase) ||
          fabsf(t.phase) > math::PI_F || t.l > 8 || t.m > t.l)
        return false;
    return true;
  }
  /** @brief Prepares the normalized field at its drift phase.
   * @param params Authored field controls.
   * @param phase Accumulated azimuthal drift, radians.
   * @return Prepared polynomial tangent field. */
  static Prepared field(const Params &params, float phase) {
    std::array<math::tangent_field::PotentialTerm, 4> terms{};
    for (size_t i = 0; i < 4; ++i) {
      const auto &t = params.terms[i];
      terms[i] = {static_cast<math::tangent_field::TermType>(t.type),
                  static_cast<math::tangent_field::Kind>(t.kind),
                  static_cast<uint8_t>(t.l),
                  static_cast<uint8_t>(t.m),
                  t.amplitude,
                  t.phase};
    }
    return Prepared(
        terms, phase,
        math::make_rotation(math::Y_AXIS, math::Vector::from_spherical(
                                              params.orientation_azimuth,
                                              params.orientation_polar)));
  }
};

/** @brief Rossby-Haurwitz zonal and travelling wave displacement. */
struct DisplaceRossbyHaurwitz
    : FieldDisplaceModel<DisplaceRossbyHaurwitz, RossbyHaurwitzParams> {
  static constexpr const char *ID = "sphere.displace.rossby-haurwitz.v1";
  static constexpr const char *NAME = "Rossby-Haurwitz Displace";
  /** @brief Checks the term topology and parameter domains.
   * @param params Authored field controls.
   * @return Whether the term configuration is admitted. */
  static bool valid_terms(const Params &params) { return params.order < 7; }
  /** @brief Prepares the normalized field at its drift phase.
   * @param params Authored field controls.
   * @param phase Accumulated azimuthal drift, radians.
   * @return Prepared polynomial tangent field. */
  static Prepared field(const Params &params, float phase) {
    using namespace math::tangent_field;
    std::array<PotentialTerm, 4> terms{};
    terms[0] = {TermType::MODE, Kind::TOROIDAL, 1, 0, -params.omega, 0};
    terms[1] = {TermType::MODE,
                Kind::TOROIDAL,
                static_cast<uint8_t>(params.order + 2),
                static_cast<uint8_t>(params.order + 1),
                params.wave * (params.order % 2 ? 1.0f : -1.0f),
                0};
    return Prepared(terms, phase);
  }
};

/** @brief Fold-aligned toroidal or poloidal invariant displacement. */
struct DisplaceSymmetric
    : FieldDisplaceModel<DisplaceSymmetric, SymmetricDisplaceParams,
                         EmptyState> {
  static constexpr const char *ID = "sphere.displace.symmetric.v1";
  static constexpr const char *NAME = "Symmetric Displace";
  /** @brief Checks the term topology and parameter domains.
   * @param params Authored field controls.
   * @return Whether the term configuration is admitted. */
  static bool valid_terms(const Params &params) {
    return params.symmetry < 9 && params.kind < 2;
  }
  /** @brief Prepares the fold-aligned invariant field.
   * @param params Authored invariant controls.
   * @return Prepared polynomial tangent field. */
  static Prepared field(const Params &params, float) {
    using namespace math::tangent_field;
    constexpr uint8_t ORDERS[] = {6, 0, 0, 0, 3, 4, 5, 6, 8};
    std::array<PotentialTerm, 4> terms{};
    const auto type = params.symmetry >= 1 && params.symmetry <= 3
                          ? static_cast<TermType>(params.symmetry)
                          : TermType::MODE;
    const uint8_t order = ORDERS[params.symmetry];
    terms[0] = {
        type, static_cast<Kind>(params.kind), order, order, params.amplitude,
        0};
    return Prepared(terms);
  }
};

} // namespace Pullback::Interp::Op
#endif
