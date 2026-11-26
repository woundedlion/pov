/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Composed effect module cases. */

#include "tests/composed_effect/context.h"

namespace hs_test::composed_effect_tests {

static_assert(
    Pullback::ComposedDetail::qualified<
        Pullback::ParameterResource<"outer_warp", Pullback::WaveShearParams,
                                    Pullback::ResourceKind::WARP>,
        CollidingWarpNames>());
static_assert(
    Pullback::ComposedDetail::qualified<
        Pullback::ParameterResource<"inner_warp", Pullback::VectorNoiseParams,
                                    Pullback::ResourceKind::WARP>,
        CollidingWarpNames>());
static_assert(CoverageFixture::Metadata::COVERAGE ==
              Pullback::ProjectionCoverageMode::WEIGHT);
static_assert(TrackedFixture::Metadata::PATH_TRACKED);
static_assert(CoverageFixture::Metadata::PATH_TRACKED);
static_assert(CoverageFixture::Metadata::OUTER_CAMERA);
static_assert(!Involution<SMALL_W, SMALL_H>::Metadata::OUTER_CAMERA);

/** @brief Direct modulation noise exposes scales beyond the baked-cube limit. */
inline void test_composed_noise_modulation_scale_limits() {
  constexpr float PROBE_SCALE = 20.0f;
  static_assert(Pattern::MAX_BAKED_NOISE_SCALE < PROBE_SCALE &&
                PROBE_SCALE <= Pattern::MAX_NOISE_SCALE);
  for (const auto evaluation :
       {Pattern::Evaluation::DIRECT, Pattern::Evaluation::BAKED_CUBE}) {
    Pullback::Modulation::Config config;
    config.evaluation = evaluation;
    config.noise.scale = PROBE_SCALE;
    const auto fields = Pullback::Modulation::control_fields(config);
    HS_EXPECT_EQ(fields[1].min, Pattern::MIN_NOISE_SCALE);
    HS_EXPECT_EQ(fields[1].max, evaluation == Pattern::Evaluation::DIRECT
                                    ? Pattern::MAX_NOISE_SCALE
                                    : Pattern::MAX_BAKED_NOISE_SCALE);
    HS_EXPECT_EQ(Pullback::Modulation::valid(config),
                 evaluation == Pattern::Evaluation::DIRECT);
  }

  reset_globals();
  AshCloud<SMALL_W, SMALL_H> effect;
  effect.init();
  const auto *scale = effect.getParameters().find("modulation.Pattern Scale");
  HS_EXPECT_TRUE(scale != nullptr);
  if (scale) {
    HS_EXPECT_EQ(scale->max, Pattern::MAX_NOISE_SCALE);
    HS_EXPECT_EQ(
        effect.updateParameter("modulation.Pattern Scale", PROBE_SCALE),
        ParamSetResult::APPLIED);
    HS_EXPECT_EQ(scale->get_requested(), PROBE_SCALE);
  }
}

/** @brief
 * KaleidoscopeFlowers mirror coordinates meet at the longitude seam across
 * presets and phases.
 */
inline void test_flowers_longitude_seam() {
  using FX = KaleidoscopeFlowers<SMALL_W, SMALL_H>;
  for (size_t preset = 0; preset < FX::PRESET_IDS.size(); ++preset) {
    const auto params = FX::preset(preset).params;
    for (int step = 0; step <= 128; ++step) {
      const auto prepared = Pullback::Warp::prepare(
          params.template get<"inner_warp">(), step / 128.0f);
      for (float latitude : {-1.4f, -0.7f, 0.0f, 0.7f, 1.4f}) {
        const auto left = Pullback::Warp::mirror_tile_coords(
            math::Complex(-math::PI_F, latitude),
            params.template get<"inner_warp">(), prepared);
        const auto right = Pullback::Warp::mirror_tile_coords(
            math::Complex(math::PI_F, latitude),
            params.template get<"inner_warp">(), prepared);
        HS_EXPECT_NEAR(left.re, right.re, 1e-5f);
        HS_EXPECT_NEAR(left.im, right.im, 1e-5f);
      }
    }
  }
}

/**
 * @brief AlienBrain holds its first preset for PRESET_DWELL_FRAMES, then
 * advances to preset 1.
 */
inline void test_alien_brain_preset_dwell() {
  using FX = AlienBrain<SMALL_W, SMALL_H>;
  reset_globals();
  FX effect;
  effect.init();

  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  for (uint16_t frame = 1; frame < FX::PRESET_DWELL_FRAMES; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(1));
}

/** @brief
 * Möbius lens angle advances at constant magnitude, pauses, and preserves
 * authored values.
 */
inline void test_mobius_grid_circular_animation() {
  using FX = MobiusGrid<SMALL_W, SMALL_H>;
  reset_globals();
  FX effect;
  effect.init();
  const auto authored = effect.serialize_parameters().params;
  effect.draw_frame();
  effect.advance_display();
  const auto first = effect.frame_for_test().resources.template get<"lens">();
  HS_EXPECT_TRUE(std::isfinite(first.magnitude));
  HS_EXPECT_GT(first.magnitude, 0.0f);
  effect.draw_frame();
  effect.advance_display();
  const auto second = effect.frame_for_test().resources.template get<"lens">();
  HS_EXPECT_EQ(second.magnitude, first.magnitude);
  HS_EXPECT_NE(first.angle, second.angle);
  verify_params_equal(effect.serialize_parameters().params, authored);
  effect.setAnimationsPaused(true);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(effect.frame_for_test().resources.template get<"lens">().angle,
               second.angle);
  effect.setAnimationsPaused(false);
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_NE(effect.frame_for_test().resources.template get<"lens">().angle,
               second.angle);
  verify_params_equal(effect.serialize_parameters().params, authored);
}

/**
 * @brief Pins exact endpoints and weighted coordinates during mapping crossfades.
 */
inline void test_palette_mapping_crossfade() {
  using namespace Pullback::Color;
  constexpr float PHASE = 0.2f;
  const auto LINEAR = PaletteMappingWeights::single(PaletteMapping::LINEAR);
  const auto CUP = PaletteMappingWeights::single(PaletteMapping::CUP);
  for (const float progress : {0.0f, 1.0f}) {
    const auto weights = PaletteMappingWeights::lerp(LINEAR, CUP, progress);
    const auto MAPPING =
        progress == 0.0f ? PaletteMapping::LINEAR : PaletteMapping::CUP;
    HS_EXPECT_EQ(weights.exact, static_cast<uint8_t>(MAPPING));
    HS_EXPECT_EQ(palette_mapping_coordinate(PHASE, weights, 1.0f, 0.0f),
                 palette_mapping_coordinate(PHASE, MAPPING, 1.0f, 0.0f));
  }
  const auto MIDPOINT = PaletteMappingWeights::lerp(LINEAR, CUP, 0.5f);
  HS_EXPECT_EQ(MIDPOINT.exact, 0xff);
  HS_EXPECT_NEAR(palette_mapping_coordinate(PHASE, MIDPOINT, 1.0f, 0.0f),
                 0.5f * PHASE + 0.5f * math::unit_cup(PHASE), 1e-6f);
  const auto BELL = PaletteMappingWeights::single(PaletteMapping::BELL);
  HS_EXPECT_NEAR(
      palette_mapping_coordinate(
          PHASE, PaletteMappingWeights::lerp(CUP, BELL, 0.5f), 1.0f, 0.0f),
      0.5f, 1e-6f);
  const auto UNCHANGED = PaletteMappingWeights::lerp(CUP, CUP, 0.3f);
  HS_EXPECT_EQ(UNCHANGED.exact, CUP.exact);
  HS_EXPECT_TRUE(UNCHANGED.values == CUP.values);
}

/** @brief Every composed effect defaults to blanking displayed columns. */
inline void test_composed_strobe_roster() {
#define HS_COMPOSED_STROBE(name)                                               \
  {                                                                            \
    HS_CONTEXT(#name);                                                         \
    reset_globals();                                                           \
    name<SMALL_W, SMALL_H> effect;                                             \
    HS_EXPECT_TRUE(effect.strobe_columns());                                   \
  }
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_STROBE)
#undef HS_COMPOSED_STROBE
}

/** @brief Millefiori pixels match its document-built operator chain. */
inline void test_millefiori_document_pipeline() {
  composed_chain_tests::verify_export<Millefiori<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
}

/** @brief AshCloud pixels match its document-built operator chain. */
inline void test_ash_cloud_document_pipeline() {
  using FX = AshCloud<SMALL_W, SMALL_H>;
  struct Frame {
    struct Temporal {
      float fraction;
    } temporal;
  };
  for (const auto &view : {math::X_AXIS, -math::Y_AXIS, math::Z_AXIS,
                           math::Vector(-0.3f, 0.4f, -0.5f).normalized()}) {
    const auto expected = lenses::dodecahedral_kaleidoscope_lens(view);
    for (float fraction : {0.0f, 1.0f}) {
      const auto actual = Pullback::Lens::DodecahedralTableKaleidoscope::apply(
          view, Frame{{fraction}});
      HS_EXPECT_EQ(actual.x, expected.x);
      HS_EXPECT_EQ(actual.y, expected.y);
      HS_EXPECT_EQ(actual.z, expected.z);
    }
  }
  composed_chain_tests::verify_export<FX>(0, DOCUMENT_CHANNEL_TOLERANCE);
}

/** @brief Promoted pixels match the authored document-built chains. */
inline void test_promoted_document_pipelines() {
  static_assert(Involution<SMALL_W, SMALL_H>::authored_instance("source") ==
                "sample");
  composed_chain_tests::verify_export<
      KaleidoscopeStainedGlass<SMALL_W, SMALL_H>>(0,
                                                  DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<KaleidoscopeHexBright<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<MobiusGrid<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<AlienCore<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<Involution<SMALL_W, SMALL_H>>(0);
  for (size_t index = 0;
       index < GridWarp<SMALL_W, SMALL_H>::authored_preset_count(); ++index)
    composed_chain_tests::verify_export<GridWarp<SMALL_W, SMALL_H>>(
        index, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<LiquidWindows<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<Iridescence<SMALL_W, SMALL_H>>(0);
  for (size_t i = 0; i < RiverStones<SMALL_W, SMALL_H>::PRESET_IDS.size(); ++i)
    composed_chain_tests::verify_export<RiverStones<SMALL_W, SMALL_H>>(i);
  for (size_t i = 0; i < Filigree<SMALL_W, SMALL_H>::PRESET_IDS.size(); ++i)
    composed_chain_tests::verify_export<Filigree<SMALL_W, SMALL_H>>(
        i, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<CosmicEyeball<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<KaleidoscopeHexOil<SMALL_W, SMALL_H>>(0);
  composed_chain_tests::verify_export<AlienBrain<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<AlienOcean<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<KaleidoscopeMandala<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<GridSpace<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<ChromaticLichen<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<MermaidSkin<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<JewelMelt<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
  composed_chain_tests::verify_export<KaleidoscopeFlowers<SMALL_W, SMALL_H>>(
      0, DOCUMENT_CHANNEL_TOLERANCE);
}
} // namespace hs_test::composed_effect_tests
