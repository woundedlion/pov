/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/operators_common.h"
#include "render/pullback/projection.h"

/**
 * @file operators_project.h
 * @brief SPHERE→PLANE crossing operator models: the projections under a
 *        wandering, spinning projection frame.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Parameter family of the projection operators.
    @details `pole-fade` is read by the attenuating projections only; the
    kernel-backed projections (peirce, bonne, airocean) saturate their value
    weight instead, so the field is inert there — a union field like the
    samples' edge-width. */
struct ProjectChainParams {
  float pole_fade = 1.0f;
  float spin_rate = 0.0f;
  float wander = 0.0f;

  static constexpr auto FIELDS = std::array{
      Field<ProjectChainParams>{"pole-fade", &ProjectChainParams::pole_fade,
                                "Pole Fade", 1.0f, 20.0f, FieldCurve::LERP},
      Field<ProjectChainParams>{
          "projection-spin-speed", &ProjectChainParams::spin_rate,
          "Projection Spin Speed", 0.0f, 0.05f, FieldCurve::LERP},
      Field<ProjectChainParams>{
          "projection-wander", &ProjectChainParams::wander, "Projection Wander",
          0.0f, 1.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<ProjectChainParams>());

/** @brief ProjectChainParams extended with the central meridian the
    meridian-consuming projections read. */
struct MeridianProjectChainParams : ProjectChainParams {
  float central_meridian = 0.0f; /**< Central meridian, in radians. */

  static constexpr auto FIELDS = concat_fields<MeridianProjectChainParams>(
      ProjectChainParams::FIELDS,
      std::array{Field<MeridianProjectChainParams>{
          "central-meridian", &MeridianProjectChainParams::central_meridian,
          "Central Meridian", 0.0f, TWO_PI_F, FieldCurve::SHORTEST_PERIODIC}});
};
static_assert(field_ids_unique<MeridianProjectChainParams>());

/** @brief Shared shape of the projection operators: the walk state, the
    frame-composed conjugate, and the per-family projection call. */
template <typename Derived, typename ParamsT> struct ProjectOpModel {
  using Input = SphereSample;
  using Output = PlaneSample;
  using Params = ParamsT;
  using State = SpatialWalkState;
  struct Prepared {
    Quaternion conjugate;
  };

  static void init(State &state, InstanceId) {
    init_walk(state, PROJECTION_WALK_SEED);
  }
  static Status migrate(State &dst, const State &src, InstanceId) {
    dst = src;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    advance_walk(state, params.wander, params.spin_rate);
  }
  static Prepared prepare(const FrameContext &ctx, const Params &,
                          const State &state) {
    return {(make_rotation(Y_AXIS, state.spin_phase) * ctx.projection_base *
             state.wander)
                .conjugate()};
  }
  static PlaneSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const Vector local = rotate(input.dir, prepared.conjugate);
    return Kernel::project(input, local, Derived::project(local, params));
  }
};

/** @brief SPHERE→PLANE crossing: the stereographic projection. */
struct ProjectStereographic
    : ProjectOpModel<ProjectStereographic, ProjectChainParams> {
  static constexpr const char *ID = "project.stereographic.v2";
  static constexpr const char *NAME = "Stereographic";

  static ProjectionResult project(const Vector &local, const Params &params) {
    return Projection::stereographic(local, params.pole_fade);
  }
};

/** @brief SPHERE→PLANE crossing: the folded sinusoidal projection. */
struct ProjectFoldedSinusoidal
    : ProjectOpModel<ProjectFoldedSinusoidal, MeridianProjectChainParams> {
  static constexpr const char *ID = "project.folded-sinusoidal.v2";
  static constexpr const char *NAME = "Folded Sinusoidal";

  static ProjectionResult project(const Vector &local, const Params &params) {
    return Projection::folded_sinusoidal(local, params.central_meridian,
                                         params.pole_fade);
  }
};

/** @brief SPHERE→PLANE crossing: the equirectangular projection. */
struct ProjectEquirectangular
    : ProjectOpModel<ProjectEquirectangular, MeridianProjectChainParams> {
  static constexpr const char *ID = "project.equirectangular.v2";
  static constexpr const char *NAME = "Equirectangular";

  static ProjectionResult project(const Vector &local, const Params &params) {
    return Projection::equirectangular(local, params.central_meridian,
                                       params.pole_fade);
  }
};

inline constexpr const char *GNOMONIC_HEMISPHERE_IDS[] = {"folded", "front",
                                                          "back"};

/** @brief Parameter family of project.gnomonic.v2. */
struct GnomonicChainParams : ProjectChainParams {
  uint8_t hemisphere =
      static_cast<uint8_t>(Projection::GnomonicHemisphere::FOLDED);

  static constexpr auto FIELDS = concat_fields<GnomonicChainParams>(
      ProjectChainParams::FIELDS, std::array<Field<GnomonicChainParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<GnomonicChainParams>{
          "hemisphere", &GnomonicChainParams::hemisphere,
          GNOMONIC_HEMISPHERE_IDS, 3,
          static_cast<uint8_t>(Projection::GnomonicHemisphere::FOLDED)},
  };
};
static_assert(field_ids_unique<GnomonicChainParams>());

/** @brief SPHERE→PLANE crossing: the gnomonic projection under a hemisphere
    topology. */
struct ProjectGnomonic : ProjectOpModel<ProjectGnomonic, GnomonicChainParams> {
  static constexpr const char *ID = "project.gnomonic.v2";
  static constexpr const char *NAME = "Gnomonic";

  static ProjectionResult project(const Vector &local, const Params &params) {
    return Projection::gnomonic(
        local, params.pole_fade,
        static_cast<Projection::GnomonicHemisphere>(params.hemisphere));
  }
};

/** Fixed kernel-projection arguments: the chain pins the square Peirce
    layout, a zero layout scroll, unit plane scale, and unconditional edge
    distances. */
inline constexpr uint8_t PEIRCE_SQUARE_LAYOUT = 1;
inline constexpr float PROJECT_COORDINATE_SCALE = 1.0f;

/** @brief SPHERE→PLANE crossing: the exact Peirce quincuncial projection on
    the square layout. */
struct ProjectPeirce
    : ProjectOpModel<ProjectPeirce, MeridianProjectChainParams> {
  static constexpr const char *ID = "project.peirce.v2";
  static constexpr const char *NAME = "Peirce";

  static ProjectionResult project(const Vector &local, const Params &params) {
    return Projection::peirce(local, params.central_meridian,
                              PEIRCE_SQUARE_LAYOUT, 0.0f, true,
                              PROJECT_COORDINATE_SCALE);
  }
};

/** @brief SPHERE→PLANE crossing: the fast square-layout Peirce
    approximation, valid only on a zero central meridian. */
struct ProjectPeirceSquareFast
    : ProjectOpModel<ProjectPeirceSquareFast, ProjectChainParams> {
  static constexpr const char *ID = "project.peirce-square-fast.v2";
  static constexpr const char *NAME = "Peirce (Fast Square)";

  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::PEIRCE_FAST_SQUARE;
  static constexpr auto METRICS = Projection::PEIRCE_FAST_SQUARE_METRICS;

  static ProjectionResult project(const Vector &local, const Params &) {
    return Projection::peirce_fast_square(local, PROJECT_COORDINATE_SCALE);
  }
};

inline constexpr const char *BONNE_HEMISPHERE_IDS[] = {"north", "south"};

/** @brief Standard parallel magnitude of the chain's Bonne projection. */
inline constexpr float BONNE_STANDARD_PARALLEL = PI_F * 0.25f;

/** @brief Parameter family of project.bonne.v2. */
struct BonneChainParams : MeridianProjectChainParams {
  uint8_t hemisphere = 0; /**< 0 north, 1 south. */

  static constexpr auto FIELDS =
      concat_fields<BonneChainParams>(MeridianProjectChainParams::FIELDS,
                                      std::array<Field<BonneChainParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<BonneChainParams>{"hemisphere",
                                      &BonneChainParams::hemisphere,
                                      BONNE_HEMISPHERE_IDS, 2, 0},
  };
};
static_assert(field_ids_unique<BonneChainParams>());

/** @brief SPHERE→PLANE crossing: the Bonne projection under a hemisphere
    topology. */
struct ProjectBonne : ProjectOpModel<ProjectBonne, BonneChainParams> {
  static constexpr const char *ID = "project.bonne.v2";
  static constexpr const char *NAME = "Bonne";

  static ProjectionResult project(const Vector &local, const Params &params) {
    const float hemisphere = params.hemisphere == 0 ? 1.0f : -1.0f;
    return Projection::bonne(local, params.central_meridian,
                             hemisphere * BONNE_STANDARD_PARALLEL,
                             PROJECT_COORDINATE_SCALE);
  }
};

/** @brief SPHERE→PLANE crossing: the airocean projection on the vertical
    layout. */
struct ProjectAirocean
    : ProjectOpModel<ProjectAirocean, MeridianProjectChainParams> {
  static constexpr const char *ID = "project.airocean.v2";
  static constexpr const char *NAME = "Airocean";

  static ProjectionResult project(const Vector &local, const Params &params) {
    return Projection::airocean(local, params.central_meridian, false, true,
                                PROJECT_COORDINATE_SCALE);
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
