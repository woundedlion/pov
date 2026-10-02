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
 * @brief SPHERE→PLANE crossing operator models with identity or spin/wander
 *        projection frames.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Projection orientation policy. */
enum class ProjectionFrame : uint8_t { IDENTITY, SPIN_WANDER };

inline constexpr const char *PROJECTION_FRAME_IDS[] = {"identity",
                                                       "spin-wander"};
static_assert(std::size(PROJECTION_FRAME_IDS) ==
              static_cast<size_t>(ProjectionFrame::SPIN_WANDER) + 1);

/** @brief Activation relation of the spin and wander rates, which the
    identity frame neither advances nor reads. */
inline constexpr TopologyGate SPIN_WANDER_FRAME_GATE{
    "frame", live_values(ProjectionFrame::SPIN_WANDER)};

/** @brief Shared projection frame topology followed by family-specific fields. */
template <typename Params, typename... Extra>
constexpr std::array<TopologyField<Params>, 1 + sizeof...(Extra)>
projection_frame_topology(const Extra &...extra) {
  return {
      TopologyField<Params>{"frame", &Params::frame, PROJECTION_FRAME_IDS,
                            static_cast<uint8_t>(ProjectionFrame::SPIN_WANDER)},
      extra...};
}

/** @brief Parameter family of the projection operators.
    @details `singularity-fade` is read by projections with a singular locus; it is
    omitted for folded sinusoidal, Bonne, and Airocean. */
struct ProjectChainParams {
  uint8_t frame = static_cast<uint8_t>(ProjectionFrame::SPIN_WANDER);
  float singularity_fade = 1.0f;
  float spin_rate = 0.0f;
  float wander = 0.0f;

  static constexpr auto FIELDS = std::array{
      Field<ProjectChainParams>{
          "singularity-fade", &ProjectChainParams::singularity_fade,
          "Singularity Fade", 1.0f, 20.0f, FieldCurve::LERP},
      Field<ProjectChainParams>{
          "projection-spin-speed", &ProjectChainParams::spin_rate,
          "Projection Spin Speed", 0.0f, 0.05f, FieldCurve::LERP,
          FieldGate::ALWAYS, SPIN_WANDER_FRAME_GATE},
      Field<ProjectChainParams>{
          "projection-wander", &ProjectChainParams::wander, "Projection Wander",
          0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
          SPIN_WANDER_FRAME_GATE},
  };
  static constexpr auto TOPOLOGY =
      projection_frame_topology<ProjectChainParams>();
};
static_assert(field_ids_unique<ProjectChainParams>());
static_assert(field_defaults_in_range<ProjectChainParams>());

/** @brief ProjectChainParams extended with the central meridian the
    meridian-consuming projections read. */
struct MeridianProjectChainParams : ProjectChainParams {
  float central_meridian = 0.0f; /**< Central meridian, in radians. */

  static constexpr auto FIELDS = concat_fields<MeridianProjectChainParams>(
      ProjectChainParams::FIELDS,
      std::array{Field<MeridianProjectChainParams>{
          "central-meridian", &MeridianProjectChainParams::central_meridian,
          "Central Meridian", 0.0f, math::TWO_PI_F,
          FieldCurve::SHORTEST_PERIODIC}});
};
static_assert(field_ids_unique<MeridianProjectChainParams>());
static_assert(appended_block_size_matches<MeridianProjectChainParams,
                                          ProjectChainParams, 4>(),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<MeridianProjectChainParams>());
static_assert(
    [] {
      constexpr MeridianProjectChainParams CHAIN_DEFAULTS{};
      constexpr Projection::ProjectionParams COMPOSED_DEFAULTS{};
      for (const auto &chain : MeridianProjectChainParams::FIELDS) {
        bool matched = false;
        for (const auto &composed : Projection::ProjectionParams::FIELDS) {
          if (std::string_view(chain.id) != composed.id)
            continue;
          matched = true;
          if (chain.min != composed.min || chain.max != composed.max ||
              chain.curve != composed.curve ||
              std::string_view(chain.name) != composed.name ||
              CHAIN_DEFAULTS.*chain.member !=
                  COMPOSED_DEFAULTS.*composed.member)
            return false;
        }
        if (!matched)
          return false;
      }
      return true;
    }(),
    "chain and composed projection parameter domains must agree");

/** @brief Projection parameters without a singularity-fade control. */
struct RegularProjectChainParams : MeridianProjectChainParams {
  static constexpr auto FIELDS = [] {
    std::array<Field<MeridianProjectChainParams>,
               MeridianProjectChainParams::FIELDS.size() - 1>
        out{};
    size_t index = 0;
    for (const auto &field : MeridianProjectChainParams::FIELDS)
      if (std::string_view(field.id) != "singularity-fade")
        out[index++] = field;
    return concat_fields<RegularProjectChainParams>(
        out, std::array<Field<RegularProjectChainParams>, 0>{});
  }();
};

static_assert(field_ids_unique<RegularProjectChainParams>());
static_assert(field_defaults_in_range<RegularProjectChainParams>());

/** @brief Prepared projection frame conjugate. */
struct ProjectOrientation {
  math::Quaternion conjugate;
};
/** @brief Projection frame with cached meridian trigonometry. */
struct MeridianProjectOrientation : ProjectOrientation {
  float meridian_cos;
  float meridian_sin;
};

/** @brief Shared projection walk state, frame conjugate, and family call. */
template <typename Derived, typename ParamsT, bool CacheMeridian = false>
struct ProjectOpModel : ValueStateModel<SpatialWalkState> {
  static constexpr bool EDGE_DISTANCE_AVAILABLE = true;
  using Input = SphereSample;
  using Output = PlaneSample;
  using Params = ParamsT;
  using Prepared = std::conditional_t<CacheMeridian, MeridianProjectOrientation,
                                      ProjectOrientation>;

  static void init(State &state, InstanceId id) {
    init_walk(state, static_cast<int32_t>(id.stable_hash));
  }
  static void validate_frame(const Params &params) {
    HS_CHECK(params.frame == static_cast<uint8_t>(ProjectionFrame::IDENTITY) ||
                 params.frame ==
                     static_cast<uint8_t>(ProjectionFrame::SPIN_WANDER),
             "projection operator: invalid frame policy");
  }
  static void advance(State &state, const Params &params) {
    validate_frame(params);
    if (params.frame == static_cast<uint8_t>(ProjectionFrame::SPIN_WANDER))
      advance_walk(state, params.wander, params.spin_rate);
  }
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    validate_frame(params);
    Prepared prepared{};
    if (params.frame == static_cast<uint8_t>(ProjectionFrame::IDENTITY))
      prepared.conjugate = math::Quaternion();
    else
      prepared.conjugate =
          (math::make_rotation(math::Y_AXIS, state.spin_phase) *
           ctx.projection_base * state.wander)
              .conjugate();
    if constexpr (CacheMeridian) {
      prepared.meridian_cos = cosf(params.central_meridian);
      prepared.meridian_sin = sinf(params.central_meridian);
    }
    return prepared;
  }
  static PlaneSample run(const SphereSample &input, const FrameContext &,
                         const Params &params, const Prepared &prepared) {
    const math::Vector local = math::rotate(input.dir, prepared.conjugate);
    if constexpr (CacheMeridian)
      return Kernel::project(input, local,
                             Derived::project(local, params, prepared));
    else
      return Kernel::project(input, local, Derived::project(local, params));
  }
};

/** @brief SPHERE→PLANE crossing: the stereographic projection. */
struct ProjectStereographic
    : ProjectOpModel<ProjectStereographic, ProjectChainParams> {
  static constexpr const char *ID = "project.stereographic.v2";
  static constexpr const char *NAME = "Stereographic";

  static ProjectionResult project(const math::Vector &local,
                                  const Params &params) {
    return Projection::stereographic(local, params.singularity_fade);
  }
};

/** @brief SPHERE→PLANE crossing: the folded sinusoidal projection. */
struct ProjectFoldedSinusoidal
    : ProjectOpModel<ProjectFoldedSinusoidal, RegularProjectChainParams> {
  static constexpr bool EDGE_DISTANCE_AVAILABLE = false;
  static constexpr const char *ID = "project.folded-sinusoidal.v2";
  static constexpr const char *NAME = "Folded Sinusoidal";

  static ProjectionResult project(const math::Vector &local,
                                  const Params &params) {
    return Projection::folded_sinusoidal(local, params.central_meridian);
  }
};

/** @brief SPHERE→PLANE crossing: the equirectangular projection. */
struct ProjectEquirectangular
    : ProjectOpModel<ProjectEquirectangular, MeridianProjectChainParams> {
  static constexpr const char *ID = "project.equirectangular.v2";
  static constexpr const char *NAME = "Equirectangular";

  static ProjectionResult project(const math::Vector &local,
                                  const Params &params) {
    return Projection::equirectangular(local, params.central_meridian,
                                       params.singularity_fade);
  }
};

inline constexpr const char *GNOMONIC_HEMISPHERE_IDS[] = {"folded", "front",
                                                          "back"};
static_assert(std::size(GNOMONIC_HEMISPHERE_IDS) ==
              static_cast<size_t>(Projection::GnomonicHemisphere::BACK) + 1);

/** @brief Parameter family of project.gnomonic.v2. */
struct GnomonicChainParams : ProjectChainParams {
  uint8_t hemisphere =
      static_cast<uint8_t>(Projection::GnomonicHemisphere::FOLDED);

  static constexpr auto TOPOLOGY =
      projection_frame_topology<GnomonicChainParams>(
          TopologyField<GnomonicChainParams>{
              "hemisphere", &GnomonicChainParams::hemisphere,
              GNOMONIC_HEMISPHERE_IDS,
              static_cast<uint8_t>(Projection::GnomonicHemisphere::FOLDED)});
};
static_assert(field_ids_unique<GnomonicChainParams>());
static_assert(
    appended_block_size_matches<GnomonicChainParams, ProjectChainParams, 1>(),
    "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<GnomonicChainParams>());

/** @brief SPHERE→PLANE crossing: the gnomonic projection under a hemisphere
    topology. */
struct ProjectGnomonic : ProjectOpModel<ProjectGnomonic, GnomonicChainParams> {
  static constexpr const char *ID = "project.gnomonic.v2";
  static constexpr const char *NAME = "Gnomonic";

  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    HS_CHECK(params.hemisphere <=
                 static_cast<uint8_t>(Projection::GnomonicHemisphere::BACK),
             "project.gnomonic: invalid hemisphere");
    return ProjectOpModel::prepare(ctx, params, state);
  }
  static ProjectionResult project(const math::Vector &local,
                                  const Params &params) {
    return Projection::gnomonic(
        local, params.singularity_fade,
        static_cast<Projection::GnomonicHemisphere>(params.hemisphere));
  }
};

/** @brief Default layout of the Peirce projection. */
inline constexpr uint8_t PEIRCE_SQUARE_LAYOUT =
    static_cast<uint8_t>(projections::PeirceLayout::SQUARE);

inline constexpr const char *BONNE_HEMISPHERE_IDS[] = {"north", "south"};
static_assert(std::size(BONNE_HEMISPHERE_IDS) == 2);

/** @brief Standard parallel magnitude of the chain's Bonne projection. */
inline constexpr float BONNE_STANDARD_PARALLEL = math::PI_F * 0.25f;

/** @brief Parameter family of project.bonne.v3. */
struct BonneChainParams : RegularProjectChainParams {
  float standard_parallel = BONNE_STANDARD_PARALLEL;
  uint8_t hemisphere = 0; /**< 0 north, 1 south. */

  static constexpr auto FIELDS = concat_fields<BonneChainParams>(
      RegularProjectChainParams::FIELDS,
      std::array{Field<BonneChainParams>{
          "standard-parallel", &BonneChainParams::standard_parallel,
          "Standard Parallel", 0.0f, math::PI_F * 0.5f, FieldCurve::LERP}});
  static constexpr auto TOPOLOGY = projection_frame_topology<BonneChainParams>(
      TopologyField<BonneChainParams>{"hemisphere",
                                      &BonneChainParams::hemisphere,
                                      BONNE_HEMISPHERE_IDS, 0});
};
static_assert(field_ids_unique<BonneChainParams>());
static_assert(appended_block_size_matches<BonneChainParams,
                                          MeridianProjectChainParams, 5>(),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<BonneChainParams>());

inline constexpr const char *AIROCEAN_LAYOUT_IDS[] = {"vertical", "horizontal"};

static_assert(std::size(AIROCEAN_LAYOUT_IDS) == 2);

/** @brief Parameter family of project.airocean.v3. */
struct AiroceanChainParams : RegularProjectChainParams {
  uint8_t layout = 0;

  static constexpr auto TOPOLOGY =
      projection_frame_topology<AiroceanChainParams>(
          TopologyField<AiroceanChainParams>{
              "layout", &AiroceanChainParams::layout, AIROCEAN_LAYOUT_IDS, 0});
};
static_assert(field_ids_unique<AiroceanChainParams>());
static_assert(appended_block_size_matches<AiroceanChainParams,
                                          RegularProjectChainParams, 1>(),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<AiroceanChainParams>());

template <typename Base> struct ScaledProjectParams : Base {
  float coordinate_scale = 1.0f;
  static constexpr auto FIELDS = concat_fields<ScaledProjectParams>(
      Base::FIELDS,
      std::array{Field<ScaledProjectParams>{
          "coordinate-scale", &ScaledProjectParams::coordinate_scale,
          "Coordinate Scale", 0.25f, 4.0f, FieldCurve::LOG_POSITIVE}});
};

inline constexpr const char *PEIRCE_LAYOUT_IDS[] = {"diamond", "square",
                                                    "horizontal", "vertical"};

static_assert(std::size(PEIRCE_LAYOUT_IDS) ==
              static_cast<size_t>(projections::PeirceLayout::VERTICAL) + 1);

struct PeirceChainParams : ScaledProjectParams<MeridianProjectChainParams> {
  float layout_scroll = 0.0f;
  uint8_t layout = PEIRCE_SQUARE_LAYOUT;
  static constexpr auto FIELDS = concat_fields<PeirceChainParams>(
      ScaledProjectParams<MeridianProjectChainParams>::FIELDS,
      std::array{Field<PeirceChainParams>{
          "layout-scroll", &PeirceChainParams::layout_scroll, "Layout Scroll",
          -1.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
          TopologyGate{"layout",
                       live_values(projections::PeirceLayout::HORIZONTAL,
                                   projections::PeirceLayout::VERTICAL)}}});
  static constexpr auto TOPOLOGY = projection_frame_topology<PeirceChainParams>(
      TopologyField<PeirceChainParams>{"layout", &PeirceChainParams::layout,
                                       PEIRCE_LAYOUT_IDS,
                                       PEIRCE_SQUARE_LAYOUT});
};

struct ProjectPeirceV3
    : ProjectOpModel<ProjectPeirceV3, PeirceChainParams, true> {
  static constexpr const char *ID = "project.peirce.v3";
  static constexpr const char *NAME = "Peirce Layout";
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    HS_CHECK(params.layout < std::size(PEIRCE_LAYOUT_IDS),
             "project.peirce: invalid layout");
    return ProjectOpModel::prepare(ctx, params, state);
  }
  static ProjectionResult project(const math::Vector &local,
                                  const Params &params,
                                  const Prepared &prepared) {
    return Projection::peirce(local, params.central_meridian, params.layout,
                              params.layout_scroll, true,
                              params.coordinate_scale, params.singularity_fade,
                              prepared.meridian_cos, prepared.meridian_sin);
  }
};

struct ProjectPeirceSquareFastV3
    : ProjectOpModel<ProjectPeirceSquareFastV3,
                     ScaledProjectParams<ProjectChainParams>> {
  static constexpr const char *ID = "project.peirce-square-fast.v3";
  static constexpr const char *NAME = "Peirce Scaled Fast Square";
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::PEIRCE_FAST_SQUARE;
  static constexpr auto METRICS = Projection::PEIRCE_FAST_SQUARE_METRICS;
  static ProjectionResult project(const math::Vector &local,
                                  const Params &params) {
    return Projection::peirce_fast_square(local, params.coordinate_scale,
                                          params.singularity_fade);
  }
};

struct ProjectBonneV3
    : ProjectOpModel<ProjectBonneV3, ScaledProjectParams<BonneChainParams>> {
  static constexpr const char *ID = "project.bonne.v3";
  static constexpr const char *NAME = "Bonne Scaled";
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    HS_CHECK(params.hemisphere < std::size(BONNE_HEMISPHERE_IDS),
             "project.bonne: invalid hemisphere");
    return ProjectOpModel::prepare(ctx, params, state);
  }
  static ProjectionResult project(const math::Vector &local,
                                  const Params &params) {
    const float hemisphere = params.hemisphere == 0 ? 1.0f : -1.0f;
    return Projection::bonne(local, params.central_meridian,
                             hemisphere * params.standard_parallel,
                             params.coordinate_scale);
  }
};

struct ProjectAiroceanV3
    : ProjectOpModel<ProjectAiroceanV3,
                     ScaledProjectParams<AiroceanChainParams>, true> {
  static constexpr const char *ID = "project.airocean.v3";
  static constexpr const char *NAME = "Airocean Scaled";
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    HS_CHECK(params.layout < std::size(AIROCEAN_LAYOUT_IDS),
             "project.airocean: invalid layout");
    return ProjectOpModel::prepare(ctx, params, state);
  }
  static ProjectionResult project(const math::Vector &local,
                                  const Params &params,
                                  const Prepared &prepared) {
    return Projection::airocean(local, params.layout == 1, true,
                                params.coordinate_scale, prepared.meridian_cos,
                                prepared.meridian_sin);
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
