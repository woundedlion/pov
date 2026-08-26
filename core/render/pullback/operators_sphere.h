/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/lens.h"
#include "render/pullback/operators_common.h"
#include "render/pullback/surface.h"

/**
 * @file operators_sphere.h
 * @brief SPHERE-endomorphism operator models: the surface displacements and
 *        the sphere-to-sphere lenses.
 */

namespace Pullback {

namespace Interp {

namespace Op {

inline constexpr const char *SURFACE_INTEGRATOR_IDS[] = {"euler", "midpoint",
                                                         "midpoint-2x"};

/** @brief Parameter family of sphere.displace.curl.v2. */
struct CurlDisplaceParams : Surface::SurfaceNoiseParams {
  uint8_t basis = static_cast<uint8_t>(::NoiseBasis::SIMPLEX);
  uint8_t integrator = static_cast<uint8_t>(Surface::Integrator::EULER);

  static constexpr auto FIELDS = concat_fields<CurlDisplaceParams>(
      Surface::SurfaceNoiseParams::FIELDS,
      std::array<Field<CurlDisplaceParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<CurlDisplaceParams>{
          "basis", &CurlDisplaceParams::basis, NOISE_BASIS_IDS, 3,
          static_cast<uint8_t>(::NoiseBasis::SIMPLEX)},
      TopologyField<CurlDisplaceParams>{
          "integrator", &CurlDisplaceParams::integrator, SURFACE_INTEGRATOR_IDS,
          3, static_cast<uint8_t>(Surface::Integrator::EULER)},
  };
};
static_assert(field_ids_unique<CurlDisplaceParams>());

/** @brief The displacement operators' prepared block: the owned noise field
    and this frame's loop point. */
struct PreparedDisplace {
  const FastNoiseLite *noise;
  Surface::PreparedLoop loop;
};

/** @brief SPHERE endomorphism: the divergence-free curl-noise displacement. */
struct DisplaceCurl : ValueStateModel<NoisePhaseState> {
  static constexpr const char *ID = "sphere.displace.curl.v2";
  static constexpr const char *NAME = "Curl Displace";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = CurlDisplaceParams;
  using Prepared = PreparedDisplace;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_noise_basis(params.basis);
    return {&state.noise, Surface::prepare(state.phase)};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &params, const Prepared &prepared) {
    return Kernel::displace(
        input,
        Surface::curl_noise(
            input.dir, *prepared.noise, static_cast<::NoiseBasis>(params.basis),
            static_cast<Surface::Integrator>(params.integrator), params.scale,
            prepared.loop.loop_offset, params.strength, true));
  }
};

/** @brief Parameter family of sphere.displace.direct.v2. */
struct DirectDisplaceParams : Surface::DirectSurfaceParams {
  uint8_t basis = static_cast<uint8_t>(::NoiseBasis::SIMPLEX);

  static constexpr auto FIELDS = concat_fields<DirectDisplaceParams>(
      Surface::DirectSurfaceParams::FIELDS,
      std::array<Field<DirectDisplaceParams>, 0>{});
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<DirectDisplaceParams>{
          "basis", &DirectDisplaceParams::basis, NOISE_BASIS_IDS, 3,
          static_cast<uint8_t>(::NoiseBasis::SIMPLEX)},
  };
};
static_assert(field_ids_unique<DirectDisplaceParams>());

/** @brief The direct displacement's prepared block: noise field, loop point
    and steering frame. */
struct PreparedDirectDisplace {
  const FastNoiseLite *noise;
  Surface::PreparedDirect direct;
};

/** @brief SPHERE endomorphism: the direction-steered noise displacement. */
struct DisplaceDirect : ValueStateModel<NoisePhaseState> {
  static constexpr const char *ID = "sphere.displace.direct.v2";
  static constexpr const char *NAME = "Direct Displace";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = DirectDisplaceParams;
  using Prepared = PreparedDirectDisplace;

  static void init(State &state, InstanceId id) { init_noise_phase(state, id); }
  static void advance(State &state, const Params &params) {
    state.phase = wrap_t(state.phase + params.speed);
  }
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
                              static_cast<::NoiseBasis>(params.basis),
                              params.scale, prepared.direct.loop_offset,
                              params.strength, prepared.direct.direction_cos,
                              prepared.direct.direction_sin, true));
  }
};

/** @brief Frame clock of sphere.displace.ripple.v2. */
struct RipplePhaseState {
  float phase = 0.0f;
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
    state.phase = fmodf(state.phase + 1.0f, params.period);
  }
  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    return Surface::prepare_ripple(params, state.phase / params.period);
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

/** @brief SPHERE endomorphism: the octant-fold glitch lens. */
struct LensGlitch : FixedLensModel<LensGlitch> {
  static constexpr const char *ID = "sphere.lens.glitch.v2";
  static constexpr const char *NAME = "Glitch Lens";
  static Vector lens(const Vector &input) { return lenses::glitch_lens(input); }
};

/** @brief SPHERE endomorphism: the longitude-dependent twist lens. */
struct LensTwist : FixedLensModel<LensTwist> {
  static constexpr const char *ID = "sphere.lens.twist.v2";
  static constexpr const char *NAME = "Twist Lens";
  static Vector lens(const Vector &input) { return lenses::twist_lens(input); }
};

/** @brief Slider range of each flat Mobius coefficient, single-sourced from
    the lens family the composed path registers over. */
inline constexpr float MOBIUS_COEFFICIENT_LIMIT =
    Lens::MobiusLensParams::COEFFICIENT_LIMIT;

/** @brief Parameter family of sphere.lens.mobius.v2: the flat coefficient
    fields the chain registers.
    @details The composed-effect path registers the same eight coefficients by
    hand over Lens::MobiusLensParams (composed_effect.h); the two spellings
    stay separate because that family nests MobiusParams and cannot carry a
    flat FIELDS table. */
struct MobiusChainParams {
  float a_re = 0.7071067811865475f;
  float a_im = 0.0f;
  float b_re = 0.0f;
  float b_im = 0.0f;
  float c_re = 0.0f;
  float c_im = 0.0f;
  float d_re = 0.7071067811865475f;
  float d_im = 0.0f;

  static constexpr auto FIELDS = std::array{
      Field<MobiusChainParams>{"mobius-a-re", &MobiusChainParams::a_re,
                               "Mobius A Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
      Field<MobiusChainParams>{"mobius-a-im", &MobiusChainParams::a_im,
                               "Mobius A Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
      Field<MobiusChainParams>{"mobius-b-re", &MobiusChainParams::b_re,
                               "Mobius B Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
      Field<MobiusChainParams>{"mobius-b-im", &MobiusChainParams::b_im,
                               "Mobius B Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
      Field<MobiusChainParams>{"mobius-c-re", &MobiusChainParams::c_re,
                               "Mobius C Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
      Field<MobiusChainParams>{"mobius-c-im", &MobiusChainParams::c_im,
                               "Mobius C Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
      Field<MobiusChainParams>{"mobius-d-re", &MobiusChainParams::d_re,
                               "Mobius D Re", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
      Field<MobiusChainParams>{"mobius-d-im", &MobiusChainParams::d_im,
                               "Mobius D Im", -MOBIUS_COEFFICIENT_LIMIT,
                               MOBIUS_COEFFICIENT_LIMIT, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<MobiusChainParams>());

/** @brief SPHERE endomorphism: the Mobius map over the flat coefficients. */
struct LensMobius : StatelessModel {
  static constexpr const char *ID = "sphere.lens.mobius.v2";
  static constexpr const char *NAME = "Mobius Lens";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = MobiusChainParams;
  struct Prepared {
    MobiusParams mobius;
  };

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &) {
    return {MobiusParams{params.a_re, params.a_im, params.b_re, params.b_im,
                         params.c_re, params.c_im, params.d_re, params.d_im}};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &, const Prepared &prepared) {
    return Kernel::lens(input, Lens::mobius(input.dir, prepared.mobius));
  }
};

enum class KaleidoscopeSymmetry : uint8_t {
  AZIMUTHAL = 0,
  TETRAHEDRAL = 1,
  OCTAHEDRAL = 2,
  DODECAHEDRAL = 3,
  TRIANGULAR_PRISM = 4,
  SQUARE_PRISM = 5,
  PENTAGONAL_PRISM = 6,
  HEXAGONAL_PRISM = 7,
  OCTAGONAL_PRISM = 8
};

inline constexpr const char *KALEIDOSCOPE_SYMMETRY_IDS[] = {
    "azimuthal",        "tetrahedral",      "octahedral",
    "dodecahedral",     "triangular-prism", "square-prism",
    "pentagonal-prism", "hexagonal-prism",  "octagonal-prism"};

/** @brief Parameter family of sphere.lens.kaleidoscope.v2: the symmetry
    topology alone. */
struct KaleidoscopeChainParams {
  uint8_t symmetry = static_cast<uint8_t>(KaleidoscopeSymmetry::AZIMUTHAL);

  static constexpr std::array<Field<KaleidoscopeChainParams>, 0> FIELDS{};
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<KaleidoscopeChainParams>{
          "symmetry", &KaleidoscopeChainParams::symmetry,
          KALEIDOSCOPE_SYMMETRY_IDS, 9,
          static_cast<uint8_t>(KaleidoscopeSymmetry::AZIMUTHAL)},
  };
};
static_assert(field_ids_unique<KaleidoscopeChainParams>());

/** @brief The symmetry switch over the shared kaleidoscope lens kernels. */
inline Vector kaleidoscope_lens(const Vector &input, uint8_t symmetry) {
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
  default:
    return lenses::kaleidoscope_lens(input);
  }
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
