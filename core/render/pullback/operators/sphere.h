/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "math/mobius.h"

#include "render/pullback/lens.h"
#include "render/pullback/operators/common.h"
#include "render/pullback/surface.h"

/**
 * @file sphere.h
 * @brief SPHERE-endomorphism operator models: the surface displacements and
 *        the sphere-to-sphere lenses.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Integrator topology values, in Surface::Integrator order. */
inline constexpr const char *SURFACE_INTEGRATOR_IDS[] = {"euler", "midpoint",
                                                         "midpoint-2x"};
static_assert(std::size(SURFACE_INTEGRATOR_IDS) ==
              static_cast<size_t>(Surface::Integrator::MIDPOINT_2X) + 1);

/** @brief Limit topology values, in math::TangentLimit order. */
inline constexpr const char *TANGENT_LIMIT_IDS[] = {"smooth", "clamp"};
static_assert(std::size(TANGENT_LIMIT_IDS) ==
              static_cast<size_t>(math::TangentLimit::CLAMP) + 1);

/** @brief Parameter family of sphere.displace.curl.v2. */
struct CurlDisplaceParams : Surface::SurfaceNoiseParams {
  /** math::NoiseBasis topology value. */
  uint8_t basis = static_cast<uint8_t>(math::NoiseBasis::SIMPLEX);
  /** Surface::Integrator topology value. */
  uint8_t integrator = static_cast<uint8_t>(Surface::Integrator::EULER);
  /** math::TangentLimit topology value. */
  uint8_t limit = static_cast<uint8_t>(math::TangentLimit::SMOOTH);

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<CurlDisplaceParams>{
          "basis", &CurlDisplaceParams::basis, NOISE_BASIS_IDS,
          static_cast<uint8_t>(math::NoiseBasis::SIMPLEX)},
      TopologyField<CurlDisplaceParams>{
          "integrator", &CurlDisplaceParams::integrator, SURFACE_INTEGRATOR_IDS,
          static_cast<uint8_t>(Surface::Integrator::EULER)},
      TopologyField<CurlDisplaceParams>{
          "limit", &CurlDisplaceParams::limit, TANGENT_LIMIT_IDS,
          static_cast<uint8_t>(math::TangentLimit::SMOOTH)},
  };
};
static_assert(field_ids_unique<CurlDisplaceParams>());
static_assert(appended_block_size_matches<CurlDisplaceParams,
                                          Surface::SurfaceNoiseParams, 3>(),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<CurlDisplaceParams>());

/** @brief The displacement operators' prepared block: the owned noise field
    and this frame's loop point. */
struct PreparedDisplace {
  const FastNoiseLite *noise; ///< The instance-owned noise field.
  Surface::PreparedLoop loop; ///< This frame's loop point.
};

/** @brief SPHERE endomorphism: the length-limited curl-noise displacement. */
struct DisplaceCurl : PhaseClockModel<NoisePhaseState> {
  static constexpr const char *ID = "sphere.displace.curl.v2";
  static constexpr const char *NAME = "Curl Displace";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = CurlDisplaceParams;
  using Prepared = PreparedDisplace;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_noise_basis(params.basis);
    HS_CHECK(params.integrator <=
                 static_cast<uint8_t>(Surface::Integrator::MIDPOINT_2X),
             "sphere.displace.curl: invalid integrator");
    HS_CHECK(params.limit <= static_cast<uint8_t>(math::TangentLimit::CLAMP),
             "sphere.displace.curl: invalid limit");
    return {&state.noise, Surface::prepare(state.phase)};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &params, const Prepared &prepared) {
    return Kernel::displace(
        input, Surface::curl_noise(
                   input.dir, *prepared.noise,
                   static_cast<math::NoiseBasis>(params.basis),
                   static_cast<Surface::Integrator>(params.integrator),
                   static_cast<math::TangentLimit>(params.limit), params.scale,
                   prepared.loop.loop_offset, params.strength, true));
  }
};

/** @brief Parameter family of sphere.displace.direct.v2. */
struct DirectDisplaceParams : Surface::DirectSurfaceParams {
  /** math::NoiseBasis topology value. */
  uint8_t basis = static_cast<uint8_t>(math::NoiseBasis::SIMPLEX);

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<DirectDisplaceParams>{
          "basis", &DirectDisplaceParams::basis, NOISE_BASIS_IDS,
          static_cast<uint8_t>(math::NoiseBasis::SIMPLEX)},
  };
};
static_assert(field_ids_unique<DirectDisplaceParams>());
static_assert(appended_block_size_matches<DirectDisplaceParams,
                                          Surface::DirectSurfaceParams, 1>(),
              "appended parameter block must have the expected rounded size");
static_assert(field_defaults_in_range<DirectDisplaceParams>());

/** @brief The direct displacement's prepared block: noise field, loop point
    and steering frame. */
struct PreparedDirectDisplace {
  const FastNoiseLite *noise;     ///< The instance-owned noise field.
  Surface::PreparedDirect direct; ///< Loop point and steering frame.
};

/** @brief SPHERE endomorphism: the direction-steered noise displacement. */
struct DisplaceDirect : PhaseClockModel<NoisePhaseState> {
  static constexpr const char *ID = "sphere.displace.direct.v2";
  static constexpr const char *NAME = "Direct Displace";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = DirectDisplaceParams;
  using Prepared = PreparedDirectDisplace;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_noise_basis(params.basis);
    return {&state.noise,
            Surface::prepare_direct(state.phase, params.direction)};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &params, const Prepared &prepared) {
    return Kernel::displace(
        input,
        Surface::direct_noise(input.dir, *prepared.noise,
                              static_cast<math::NoiseBasis>(params.basis),
                              params.scale, prepared.direct.loop_offset,
                              params.strength, prepared.direct.direction_cos,
                              prepared.direct.direction_sin, true));
  }
};

/** @brief Frame clock of sphere.displace.ripple.v2. */
struct RipplePhaseState {
  float phase = 0.0f; ///< Frames into the current period.
};

/** @brief SPHERE endomorphism: a periodically expanding Ricker-wave ripple. */
struct DisplaceRipple : ValueStateModel<RipplePhaseState> {
  static constexpr const char *ID = "sphere.displace.ripple.v2";
  static constexpr const char *NAME = "Ripple Displace";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = Surface::PeriodicRippleParams;
  using Prepared = Surface::PreparedRipple;

  static void advance(State &state, const Params &params) {
    Surface::advance_ripple_phase(state.phase, params);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Surface::prepare_ripple(params,
                                   Surface::ripple_cycle(state.phase, params));
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &, const Prepared &prepared) {
    return Kernel::displace(
        input, Surface::periodic_ripple(input.dir, prepared, true));
  }
};

/** @brief Shared shape of the parameterless lens operators. */
template <typename Derived> struct FixedLensModel : StatelessModel {
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = Lens::NoLensParams;
  struct Prepared {};

  static Prepared prepare(const FrameContext &, const Params &,
                          const StatelessModel::State &) {
    return {};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &, const Prepared &) {
    return Kernel::lens(input, Derived::lens(input.dir));
  }
};

/**
 * @brief SPHERE endomorphism: the glitch lens (latitude doubling, azimuth
 *        tripling).
 */
struct LensGlitch : FixedLensModel<LensGlitch> {
  static constexpr const char *ID = "sphere.lens.glitch.v2";
  static constexpr const char *NAME = "Glitch Lens";
  /**
   * @brief The glitch lens map.
   * @param input Unit direction.
   * @return The lensed direction.
   */
  static math::Vector lens(const math::Vector &input) {
    return lenses::glitch_lens(input);
  }
};

/** @brief Parameter family of sphere.lens.twist.v2. */
struct TwistChainParams {
  float twist_rate = lenses::TWIST_RATE; /**< Radians per unit height. */

  static constexpr auto FIELDS = std::array{
      Field<TwistChainParams>{"twist-rate", &TwistChainParams::twist_rate,
                              "Twist Rate", -12.0f, 12.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<TwistChainParams>());
static_assert(field_defaults_in_range<TwistChainParams>());

/** @brief SPHERE endomorphism: the height-dependent twist lens. */
struct LensTwist : StatelessModel {
  static constexpr const char *ID = "sphere.lens.twist.v2";
  static constexpr const char *NAME = "Twist Lens";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = TwistChainParams;
  struct Prepared {};

  static Prepared prepare(const FrameContext &, const Params &, const State &) {
    return {};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &params, const Prepared &) {
    return Kernel::lens(input,
                        lenses::twist_lens(input.dir, params.twist_rate));
  }
};

/** @brief Slider range of each flat Mobius coefficient. */
inline constexpr float MOBIUS_COEFFICIENT_LIMIT =
    Lens::MobiusLensParams::COEFFICIENT_LIMIT;

/** @brief Parameter family of sphere.lens.mobius.v2: the flat coefficient
    fields the chain registers. */
struct MobiusChainParams {
  float a_re = 0.7071067811865475f; /**< Re a of (az + b) / (cz + d). */
  float a_im = 0.0f;                /**< Im a. */
  float b_re = 0.0f;                /**< Re b. */
  float b_im = 0.0f;                /**< Im b. */
  float c_re = 0.0f;                /**< Re c. */
  float c_im = 0.0f;                /**< Im c. */
  float d_re = 0.7071067811865475f; /**< Re d. */
  float d_im = 0.0f;                /**< Im d. */

  static constexpr auto FIELDS = std::array{
      Field<MobiusChainParams>{"mobius-a-re", &MobiusChainParams::a_re,
                               "Mobius A Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
      Field<MobiusChainParams>{"mobius-a-im", &MobiusChainParams::a_im,
                               "Mobius A Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
      Field<MobiusChainParams>{"mobius-b-re", &MobiusChainParams::b_re,
                               "Mobius B Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
      Field<MobiusChainParams>{"mobius-b-im", &MobiusChainParams::b_im,
                               "Mobius B Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
      Field<MobiusChainParams>{"mobius-c-re", &MobiusChainParams::c_re,
                               "Mobius C Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
      Field<MobiusChainParams>{"mobius-c-im", &MobiusChainParams::c_im,
                               "Mobius C Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
      Field<MobiusChainParams>{"mobius-d-re", &MobiusChainParams::d_re,
                               "Mobius D Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
      Field<MobiusChainParams>{"mobius-d-im", &MobiusChainParams::d_im,
                               "Mobius D Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::SNAP},
  };
};
static_assert(field_ids_unique<MobiusChainParams>());
static_assert(field_defaults_in_range<MobiusChainParams>());

/** @brief SPHERE endomorphism: the Mobius map over the flat coefficients. */
struct LensMobius : StatelessModel {
  static constexpr const char *ID = "sphere.lens.mobius.v2";
  static constexpr const char *NAME = "Mobius Lens";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = MobiusChainParams;
  /** @brief The frame's validated coefficients. */
  struct Prepared {
    math::MobiusParams mobius; ///< Nondegenerate Mobius coefficients.
  };

  /**
   * @brief Gathers the flat coefficient fields.
   * @param params Param block.
   * @return The Mobius coefficients.
   */
  static math::MobiusParams coefficients(const Params &params) {
    return {params.a_re, params.a_im, params.b_re, params.b_im,
            params.c_re, params.c_im, params.d_re, params.d_im};
  }

  /**
   * @brief Rejects degenerate coefficients.
   * @param params Param block.
   * @return The degenerate-map warning, or null.
   */
  static const char *validate(const Params &params) {
    const math::MobiusParams mobius = coefficients(params);
    return Lens::MobiusLensParams::nondegenerate(mobius)
               ? nullptr
               : Lens::MobiusLensParams::DEGENERATE_WARNING;
  }

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    const math::MobiusParams mobius = coefficients(params);
    HS_CHECK(Lens::MobiusLensParams::nondegenerate(mobius),
             "sphere.lens.mobius: degenerate coefficients");
    return {mobius};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &, const Prepared &prepared) {
    return Kernel::lens(input, Lens::mobius(input.dir, prepared.mobius));
  }
};

/** @brief Symmetry group of the kaleidoscope lens. */
enum class KaleidoscopeSymmetry : uint8_t {
  AZIMUTHAL = 0,        ///< Azimuthal wedge fold.
  TETRAHEDRAL = 1,      ///< Tetrahedral group.
  OCTAHEDRAL = 2,       ///< Octahedral group.
  DODECAHEDRAL = 3,     ///< Icosahedral/dodecahedral group.
  TRIANGULAR_PRISM = 4, ///< 3-fold prism group.
  SQUARE_PRISM = 5,     ///< 4-fold prism group.
  PENTAGONAL_PRISM = 6, ///< 5-fold prism group.
  HEXAGONAL_PRISM = 7,  ///< 6-fold prism group.
  OCTAGONAL_PRISM = 8   ///< 8-fold prism group.
};

/** @brief Symmetry topology values, in KaleidoscopeSymmetry order. */
inline constexpr const char *KALEIDOSCOPE_SYMMETRY_IDS[] = {
    "azimuthal",        "tetrahedral",      "octahedral",
    "dodecahedral",     "triangular-prism", "square-prism",
    "pentagonal-prism", "hexagonal-prism",  "octagonal-prism"};
static_assert(std::size(KALEIDOSCOPE_SYMMETRY_IDS) ==
              static_cast<size_t>(KaleidoscopeSymmetry::OCTAGONAL_PRISM) + 1);

/** @brief Parameter family of sphere.lens.kaleidoscope.v2: the symmetry
    topology alone. */
struct KaleidoscopeChainParams {
  /** KaleidoscopeSymmetry topology value. */
  uint8_t symmetry = static_cast<uint8_t>(KaleidoscopeSymmetry::AZIMUTHAL);

  static constexpr std::array<Field<KaleidoscopeChainParams>, 0> FIELDS{};
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<KaleidoscopeChainParams>{
          "symmetry", &KaleidoscopeChainParams::symmetry,
          KALEIDOSCOPE_SYMMETRY_IDS,
          static_cast<uint8_t>(KaleidoscopeSymmetry::AZIMUTHAL)},
  };
};
static_assert(field_ids_unique<KaleidoscopeChainParams>());
static_assert(field_defaults_in_range<KaleidoscopeChainParams>());

/** @brief The symmetry switch over the shared kaleidoscope lens kernels. */
inline math::Vector kaleidoscope_lens(const math::Vector &input,
                                      uint8_t symmetry) {
  switch (static_cast<KaleidoscopeSymmetry>(symmetry)) {
  case KaleidoscopeSymmetry::TETRAHEDRAL:
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::TETRAHEDRAL_MIRRORS);
  case KaleidoscopeSymmetry::OCTAHEDRAL:
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::OCTAHEDRAL_MIRRORS);
  case KaleidoscopeSymmetry::DODECAHEDRAL:
    return lenses::dodecahedral_kaleidoscope_lens(input);
  case KaleidoscopeSymmetry::TRIANGULAR_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::TRIANGULAR_PRISM_MIRRORS);
  case KaleidoscopeSymmetry::SQUARE_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::SQUARE_PRISM_MIRRORS);
  case KaleidoscopeSymmetry::PENTAGONAL_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::PENTAGONAL_PRISM_MIRRORS);
  case KaleidoscopeSymmetry::HEXAGONAL_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::HEXAGONAL_PRISM_MIRRORS);
  case KaleidoscopeSymmetry::OCTAGONAL_PRISM:
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::OCTAGONAL_PRISM_MIRRORS);
  case KaleidoscopeSymmetry::AZIMUTHAL:
    break;
  }
  return lenses::kaleidoscope_lens(input);
}

/** @brief SPHERE endomorphism: the kaleidoscope lens under a symmetry
    topology. */
struct LensKaleidoscope : StatelessModel {
  static constexpr const char *ID = "sphere.lens.kaleidoscope.v2";
  static constexpr const char *NAME = "Kaleidoscope Lens";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = KaleidoscopeChainParams;
  struct Prepared {};

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    HS_CHECK(params.symmetry <=
                 static_cast<uint8_t>(KaleidoscopeSymmetry::OCTAGONAL_PRISM),
             "sphere.lens.kaleidoscope: invalid symmetry");
    return {};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &params, const Prepared &) {
    return Kernel::lens(input, kaleidoscope_lens(input.dir, params.symmetry));
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
