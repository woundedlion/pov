/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file AshCloud.h
 * @brief Composed pullback effect:
 *        a softly cut lattice curled across dodecahedral facets.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Canonical dodecahedral lens with validated interior matrix lookup. */
struct AshCloudDodecahedralLens : Pullback::Lens::DodecahedralKaleidoscope {
  /** @brief Preserves literal endpoints and validates interior matrix folds.
   * @tparam FrameState Prepared column context.
   * @param input Unit view direction.
   * @param frame Column context.
   * @return Canonical chamber direction. */
  template <typename FrameState>
  HS_HOT_INLINE static math::Vector apply(const math::Vector &input,
                                          const FrameState &frame) {
    if constexpr (requires { frame.temporal.fraction; })
      if (frame.temporal.fraction > 0.0f && frame.temporal.fraction < 1.0f) {
        const auto *matrices = lenses::DodecahedralTable::GROUP.values.data();
        const uint8_t *faces = nullptr;
        if constexpr (requires { frame.dodecahedral_matrices; })
          if (frame.dodecahedral_matrices)
            matrices = frame.dodecahedral_matrices;
        if constexpr (requires { frame.dodecahedral_faces; })
          faces = frame.dodecahedral_faces;
        return lenses::dodecahedral_table_lens(input, matrices, faces);
      }
    return lenses::dodecahedral_kaleidoscope_lens(input);
  }
};

/** @brief Euler integration with an inline qualified gradient core. */
struct AshCloudInlineEuler : Pullback::Surface::Euler {
  static constexpr bool INLINE_UNIT_GRADIENT =
      true; ///< Inlines unit-frequency gradients.
};

/** @brief Pullback stage spec for `AshCloud`. */
struct AshCloudSpec : Pullback::Spec {
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr bool CACHE_DODECAHEDRAL_FACES =
      true; ///< Copies cube classifiers into RAM.
  static constexpr bool CACHE_DODECAHEDRAL_MATRICES =
      true; ///< Copies reflection matrices into RAM.
  static constexpr bool TEMPORAL_COLOR_PHASE = true; ///< Shared color tail.
  /// Caches interval palette blends in OKLab.
  static constexpr bool PALETTE_LAB_CACHE = true;
  /// Caches final palette colors by value and hue.
  static constexpr bool PALETTE_RGB_CACHE = true;
  static constexpr float PALETTE_RGB_CHROMA =
      .314f;                                      ///< Qualified source chroma.
  static constexpr bool PHASE_NOISE_CACHE = true; ///< Cached chamber noise.
  static constexpr bool PALETTE_RGB_MIDPOINT_ONLY =
      true; ///< Uses one admitted RGB midpoint.
  static constexpr bool CACHE_FULL_GAMUT_LUT =
      true; ///< Persistent gamut lookup.
  static constexpr bool RIPPLE_RADIAL_TABLE =
      true; ///< Cached ripple transfers.
  static constexpr bool PERSISTENT_PREPARED_FRAME =
      true; ///< Arena-owned frame.
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Stage chain: camera rotation, burst ripples, surface curl,
   * dodecahedral kaleidoscope lens, projection, lattice source, value cutout,
   * generated palette, noise and displacement hue modulation.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Displace<
          Pullback::Surface::BurstRipple<Pullback::SurfaceProvider<
              B, Pullback::SphereAnimation::BurstRippleParams, true,
              "displace">>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Displace<Pullback::Surface::CurlNoise<
              Pullback::SurfaceProvider<B, Pullback::SurfaceNoiseParams, true>,
              math::NoiseBasis::SIMPLEX, AshCloudInlineEuler,
              math::TangentLimit::SMOOTH>>,
          Pullback::Stage::Lens<AshCloudDodecahedralLens>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::STEREOGRAPHIC, B>::Type>>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::LatticeSourceParams, B>::Type,
                              Pullback::Weight::Projection,
                              typename Pullback::CoveragePolicyFor<
                                  Pullback::ProjectionCoverageMode::WEIGHT, B,
                                  Pullback::CutoutValueParams>::Type>,
      Pullback::Stage::ApplyCoverage<Pullback::ValueCoverage::ValueCutout<
          Pullback::ValueProvider<B, Pullback::CutoutValueParams>>>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, Pullback::Color::BrightnessEnvelope::NONE,
                                  PALETTE_LAB_CACHE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::DIRECT, TEMPORAL_COLOR_PHASE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::DISPLACEMENT,
              Pattern::Evaluation::DIRECT, TEMPORAL_COLOR_PHASE>>>>;
};
/// Parameter block derived from `AshCloudSpec`.
using AshCloudParams = Pullback::ParamsFor<AshCloudSpec>;

/**
 * @brief A softly cut lattice curled across dodecahedral facets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AshCloud
    : public Pullback::ComposedEffect<W, H, AshCloud<W, H>, AshCloudSpec> {

  using Base = Pullback::ComposedEffect<W, H, AshCloud<W, H>, AshCloudSpec>;

public:
  /// Parameter block derived from `AshCloudSpec`.
  using Params = AshCloudParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::AshCloud::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "40605c23a59db06f66579a3b652c5ba45ab04e8f54c1e43f05e7453cf9bb80ff";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "38bddb9ffc08b44f637e3af83ecb9daaeafa83958523e4f94c8abc442c291e62";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "ash-cloud"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  /// Outer camera spin about +Y, in radians per frame.
  static constexpr float CAMERA_SPIN_RATE = 0.01975f;
  // clang-format on
  // End generated identity.
  /// Snapshot schema version; changes with the `Params` layout.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 5;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().mapping_frequency = 3.666f;
    value.template get<"color">().mapping_phase = 0.598f;
    value.template get<"color">().palette_chroma = 0.314f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.856f;
    value.template get<"color">().phase_oscillation_speed = 0.00448f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"value">().cutout_softness = 0.31635937f;
    value.template get<"value">().cutout_threshold = 0.365f;
    value.template get<"displace">().burst_rate = 0.015625f;
    value.template get<"displace">().center_azimuth = 0.0f;
    value.template get<"displace">().center_polar = 1.5707964f;
    value.template get<"displace">().decay = 0.1f;
    value.template get<"displace">().period = 80.0f;
    value.template get<"displace">().strength = 0.15f;
    value.template get<"displace">().thickness = 0.7f;
    value.template get<"modulate">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"modulate">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"modulate">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"modulate">().modulation.grid.axis.x = 0.0f;
    value.template get<"modulate">().modulation.grid.axis.y = 1.0f;
    value.template get<"modulate">().modulation.grid.axis.z = 0.0f;
    value.template get<"modulate">().modulation.grid.frequency = 1.0f;
    value.template get<"modulate">().modulation.grid.frequency_b = 1.0f;
    value.template get<"modulate">().modulation.grid.offset = 0.0f;
    value.template get<"modulate">().modulation.grid.offset_b = 0.0f;
    value.template get<"modulate">().modulation.grid.speed = 0.0f;
    value.template get<"modulate">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"modulate">().modulation.interference.axis_b.y = 0.0f;
    value.template get<"modulate">().modulation.interference.axis_b.z = 0.0f;
    value.template get<"modulate">().modulation.interference.axis.x = 0.0f;
    value.template get<"modulate">().modulation.interference.axis.y = 1.0f;
    value.template get<"modulate">().modulation.interference.axis.z = 0.0f;
    value.template get<"modulate">().modulation.interference.frequency = 1.0f;
    value.template get<"modulate">().modulation.interference.frequency_b = 1.0f;
    value.template get<"modulate">().modulation.interference.offset = 0.0f;
    value.template get<"modulate">().modulation.interference.offset_b = 0.0f;
    value.template get<"modulate">().modulation.interference.speed = 0.0f;
    value.template get<"modulate">().modulation.mechanic.amount = -1.3882353f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
    value.template get<"modulate">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate">().modulation.noise.scale = 1.0f;
    value.template get<"modulate">().modulation.noise.speed = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate">().modulation.sine.offset = 0.0f;
    value.template get<"modulate">().modulation.sine.speed = 0.0f;
    value.template get<"modulation">().modulation.mechanic.amount = -3.0823529f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulation">().modulation.noise.scale = 2.4029531f;
    value.template get<"modulation">().modulation.noise.speed = 0.00015f;
    value.template get<"projection">().spin_rate = 0.00175f;
    value.template get<"projection">().wander = 0.009f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().lattice_cell_scale = 0.98971874f;
    value.template get<"source">().lattice_radius = 0.37757313f;
    value.template get<"source">().lattice_shape_blend = 0.0f;
    value.template get<"source">().lattice_softness = 0.34563965f;
    value.template get<"surface">().scale = 5.7742186f;
    value.template get<"surface">().speed = 0.00040625f;
    value.template get<"surface">().strength = 0.043f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Evaluates the canonical composed pixel pipeline.
   * @param view Unit view direction.
   * @param frame Prepared column state.
   * @return Straight-alpha pixel color. */
  HS_HOT_FLASH_MEMBER static Color4 shade(const math::Vector &view,
                                          const typename Base::Frame &frame) {
    return Base::RenderPipeline::evaluate(view, frame.ctx, frame.prepared);
  }

  /** @brief Shades a column admitted to the shared cached color path.
   * @param rotated Camera-rotated unit direction.
   * @param frame Prepared eligible column.
   * @return Straight-alpha pixel color. */
  HS_HOT_ITCM_MEMBER static Color4
  shade_shared_rotated(const math::Vector &rotated,
                       const typename Base::Frame &frame) {
    using Pipeline = typename Base::RenderPipeline;
    static_assert(Pipeline::NODE_COUNT == 8);
    const auto RIPPLED = Pipeline::template NodeAt<1>::run(
        Pullback::SphereSample{rotated, 0.0f}, frame.ctx,
        std::get<1>(frame.prepared));
    const auto PROJECTED = Pipeline::template NodeAt<2>::run(
        RIPPLED, frame.ctx, std::get<2>(frame.prepared));
    const auto SAMPLED = Pipeline::template NodeAt<3>::run(
        PROJECTED, frame.ctx, std::get<3>(frame.prepared));
    const auto COVERED = Pipeline::template NodeAt<4>::run(
        SAMPLED, frame.ctx, std::get<4>(frame.prepared));
    const auto &BASE = std::get<5>(frame.prepared);
    const auto &STATE = BASE.current;
    const float VALUE = Pullback::Color::palette_mapping_coordinate(
        COVERED.value, STATE.mapping, STATE.mapping_frequency,
        STATE.mapping_offset);
    const float OPACITY =
        COVERED.coverage *
        hs::lerp(STATE.opacity_low, STATE.opacity_high, COVERED.value);
    float hue = 0.0f;
    const auto apply_hue = [&]<size_t Index>() HS_HOT_INLINE {
      const auto &PREPARED = std::get<Index>(frame.prepared);
      if (!PREPARED.enabled || PREPARED.mechanic.amount == 0.0f)
        return;
      using Policy =
          typename Pipeline::template NodeAt<Index>::Descriptor::ColorPolicy;
      const float SCALAR = Policy::sample_phase_scalar(COVERED, PREPARED);
      hue = math::wrap_t(hue + math::wrap_t(PREPARED.mechanic.amount * SCALAR));
    };
    apply_hue.template operator()<6>();
    apply_hue.template operator()<7>();
    return Color4(PaletteRgbHueLut::sample(BASE.lab, VALUE, hue), OPACITY);
  }

  /** @brief Column shader with a camera-rotated meridian basis. */
  struct ColumnShader {
    const typename Base::Frame *frame; ///< Borrowed prepared column.
    math::Vector equator;              ///< Rotated equatorial direction.
    math::Vector up;                   ///< Rotated polar direction.
    float cosine;                      ///< Column azimuth cosine.
    float sine;                        ///< Column azimuth sine.
    bool shared;   ///< Whether the column admits the shared cache.
    bool factored; ///< Whether the physical camera grid was qualified.

    /** @brief Shades one latitude with its selected color path.
     * @param row_sine Sine of the row polar angle.
     * @param row_cosine Cosine of the row polar angle.
     * @return Straight-alpha color. */
    Color4 operator()(float row_sine, float row_cosine) const {
      const math::Vector VIEW(row_sine * cosine, row_cosine, row_sine * sine);
      if (!shared)
        return AshCloud::shade(VIEW, *frame);
      math::Vector rotated;
      if (factored)
        rotated = row_sine * equator + row_cosine * up;
      else
        rotated = Base::RenderPipeline::template NodeAt<0>::run(
                      Pullback::SphereSample{VIEW, 0.0f}, frame->ctx,
                      std::get<0>(frame->prepared))
                      .dir;
      return AshCloud::shade_shared_rotated(rotated, *frame);
    }
  };

  /** @brief Whether the current grid admits factored camera arithmetic.
   * @return True for the qualified physical display grid. */
  static bool camera_column_admitted() {
    return W == 288 && H == 144 &&
           math::DISPLAY_NORTH_PHI == .02f * math::PI_F &&
           math::DISPLAY_SOUTH_PHI == .98f * math::PI_F &&
           math::DisplayGeometry<H>::NORTH_PHI == .02f * math::PI_F &&
           math::DisplayGeometry<H>::SOUTH_PHI == .98f * math::PI_F;
  }

  /** @brief Selects the shared shader and factors its camera once per column.
   * @param frame Prepared column.
   * @param x Absolute column index.
   * @return Shader borrowing this column's frame. */
  static ColumnShader prepare_column_shader(const typename Base::Frame &frame,
                                            int x) {
    using Pipeline = typename Base::RenderPipeline;
    const auto &BASE = std::get<5>(frame.prepared);
    const auto hue_only = [](const auto &prepared) {
      return !prepared.enabled || prepared.mechanic.amount == 0.0f ||
             prepared.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION;
    };
    const bool SHARED =
        BASE.current.brightness_envelope ==
            Pullback::Color::BrightnessEnvelope::NONE &&
        hue_only(std::get<6>(frame.prepared)) &&
        hue_only(std::get<7>(frame.prepared)) && BASE.lab.low &&
        BASE.lab.rgb_low && BASE.lab.rgb_high && BASE.fraction > 0.0f &&
        BASE.fraction < 1.0f &&
        Pipeline::template NodeAt<6>::phase_ready(
            std::get<6>(frame.prepared)) &&
        Pipeline::template NodeAt<7>::phase_ready(std::get<7>(frame.prepared));
    ColumnShader shader{&frame,
                        {},
                        {},
                        math::TrigLUT<W, H>::cos_theta(x),
                        math::TrigLUT<W, H>::sin_theta[x],
                        SHARED,
                        SHARED && camera_column_admitted()};
    if (shader.factored) {
      shader.equator =
          Pipeline::template NodeAt<0>::run(
              Pullback::SphereSample{{shader.cosine, 0.0f, shader.sine}, 0.0f},
              frame.ctx, std::get<0>(frame.prepared))
              .dir;
      shader.up = Pipeline::template NodeAt<0>::run(
                      Pullback::SphereSample{{0.0f, 1.0f, 0.0f}, 0.0f},
                      frame.ctx, std::get<0>(frame.prepared))
                      .dir;
    }
    return shader;
  }

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, AshCloud<W, H>,
                                    AshCloudSpec>::choreography_hooks()
        .template extend<
            AshCloud,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &AshCloud::initial_params>>();
  }
};
