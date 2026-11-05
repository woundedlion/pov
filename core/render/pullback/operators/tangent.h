/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER
#include "core/math/tangent_field.h"
#include "core/render/pullback/tangent_surface.h"
#include "core/render/pullback/operators/sphere.h"

/** @file tangent.h
 * @brief Fixed-step polynomial tangent-field operators. */

namespace Pullback::Interp::Op {

/** @brief Scalar-potential topology names. */
using Surface::POTENTIAL_TERM_IDS;
/** @brief Helmholtz component topology names. */
using Surface::POTENTIAL_KIND_IDS;
/** @brief Harmonic degree and order values, indexed by the integer value. */
using Surface::HARMONIC_ORDER_IDS;
/** @brief Rossby wave-order values; the stored index represents R = index + 1. */
inline constexpr auto &ROSSBY_ORDER_IDS = Surface::AUTHORED_ROSSBY_ORDER_IDS;
/** @brief Rejected field configuration diagnostic. */
inline constexpr const char *FIELD_RANGE_WARNING =
    "Field parameters must be finite and within their scalar and harmonic domains.";

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

/** @brief Canonical scalar-potential term controls. */
using PotentialTermParams = Surface::PotentialTermParams;
/** @brief Canonical four-term polynomial field controls. */
using PotentialDisplaceParams = Surface::PotentialParams;

/** @brief Rossby-Haurwitz field controls. */
using RossbyHaurwitzParams = Surface::RossbyHaurwitzParams;

/** @brief Rotating invariant field controls. */
using SymmetricDisplaceParams = Surface::SymmetricParams;

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
  /// Interpolation admission contract for field parameter domains.
  static constexpr const char *ADMISSIBILITY_CONVEXITY =
      "Scalar fields retain their declared domains; harmonic topology snaps.";
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
  /** @brief Admits field families without nested potential terms.
   * @return True. */
  static constexpr bool valid_terms(const Params &) { return true; }
  /** @brief Checks field scalar and harmonic domains.
   * @param params Authored field controls.
   * @return Null when admitted, otherwise the range diagnostic. */
  static const char *validate(const Params &params) {
    if (!Fields::valid(params) || !Derived::valid_terms(params))
      return FIELD_RANGE_WARNING;
    return nullptr;
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
    return Surface::step_tangent_field(p, field, strength, integrator, true);
  }
};

/** @brief Four scalar-potential terms assembled into one displacement. */
struct DisplacePotential
    : FieldDisplaceModel<DisplacePotential, PotentialDisplaceParams> {
  static constexpr const char *ID = "sphere.displace.potential.v1";
  static constexpr const char *NAME = "Potential Lens";
  /** @brief Checks the term topology and parameter domains.
   * @param params Authored field controls.
   * @return Whether the term configuration is admitted. */
  static bool valid_terms(const Params &params) {
    return Surface::valid_terms(params);
  }
  /** @brief Prepares the normalized field at its drift phase.
   * @param params Authored field controls.
   * @param phase Accumulated azimuthal drift, radians.
   * @return Prepared polynomial tangent field. */
  static Prepared field(const Params &params, float phase) {
    return Surface::potential_field(params, phase);
  }
};

/** @brief Rossby-Haurwitz zonal and travelling wave displacement. */
struct DisplaceRossbyHaurwitz
    : FieldDisplaceModel<DisplaceRossbyHaurwitz, RossbyHaurwitzParams> {
  static constexpr const char *ID = "sphere.displace.rossby-haurwitz.v1";
  static constexpr const char *NAME = "Rossby-Haurwitz Lens";
  /** @brief Prepares the normalized field at its drift phase.
   * @param params Authored field controls.
   * @param phase Accumulated azimuthal drift, radians.
   * @return Prepared polynomial tangent field. */
  static Prepared field(const Params &params, float phase) {
    return Surface::rossby_haurwitz_field(params.omega, params.wave,
                                          params.order, phase);
  }
};

/** @brief Rotating toroidal or poloidal invariant displacement. */
struct DisplaceSymmetric
    : FieldDisplaceModel<DisplaceSymmetric, SymmetricDisplaceParams> {
  static constexpr const char *ID = "sphere.displace.symmetric.v2";
  static constexpr const char *NAME = "Symmetric Lens";
  static void advance(State &state, const Params &params) {
    state.phase = fmodf(state.phase + params.speed, math::TWO_PI_F);
  }
  /** @brief Prepares the rotated invariant field.
   * @param params Authored invariant controls.
   * @param phase Accumulated frame rotation, radians.
   * @return Prepared polynomial tangent field. */
  static Prepared field(const Params &params, float phase) {
    return Surface::symmetric_field(params.symmetry, params.kind, phase);
  }
};

} // namespace Pullback::Interp::Op
#endif
