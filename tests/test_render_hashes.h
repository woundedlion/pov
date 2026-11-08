/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Render hashes of every roster effect at every preset, for
 *        tools/render_diff.sh.
 */
#pragma once

#include "core/render/canvas.h"
#include "targets/wasm/effect_factory.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cstdio>
#include <memory>

namespace hs_test {
namespace render_hashes_tests {

/** @brief Render width in pixels. */
constexpr int HASH_W = 96;
/** @brief Render height in pixels. */
constexpr int HASH_H = 20;
/** @brief Frames folded into each preset's hash. */
constexpr int HASH_FRAMES = 64;

/**
 * @brief Renders one roster effect at one preset with daydream's load
 *        settings and prints its hash when HS_PRINT_RENDER_HASHES is set.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param entry Factory entry to build.
 * @param preset Preset selected before the first frame; ignored when the
 *        effect authors no presets.
 * @details Seeds the effect stream from the entry's identity, keeps default
 *          parameters, renders on one frame timestamp (temporal sampling off)
 *          under a pinned frame clock, and requires a lit frame.
 */
template <int W, int H>
inline void check_preset_render(const FactoryEntry &entry, size_t preset) {
  HS_CONTEXT(entry.identity.display_name.data());
  reset_globals();
  Spacetime::set_flat_surface(true);
  hs::seed_effect_stream(entry.identity);
  pin_frame_clock(0);
  std::unique_ptr<Effect> effect = entry.creator();
  HS_EXPECT_TRUE(effect != nullptr);
  if (!effect)
    return;
  effect->init();
  HS_EXPECT_EQ(effect->getPresetCount(), entry.preset_count);
  if (entry.preset_count > 0)
    HS_EXPECT_TRUE(effect->selectPreset(preset));

  uint64_t hash = HASH_SEED64;
  bool lit = false;
  for (int f = 0; f < HASH_FRAMES; ++f) {
    pin_frame_clock(f);
    effect->draw_frame();
    effect->advance_display();
    hash = frame_hash<W, H>(*effect, hash);
    lit = lit || count_lit_region<W, H>(*effect) > 0;
  }
  HS_EXPECT_TRUE(lit);

  char scenario[24];
  std::snprintf(scenario, sizeof(scenario), "p%zu", preset);
  print_render_hash(entry.identity.display_name.data(), scenario, hash);
  effect.reset();
  hs::clear_mock_time();
}

/**
 * @brief Hashes every factory entry at every preset when
 *        HS_PRINT_RENDER_HASHES is set; otherwise checks only the roster size.
 */
inline void test_roster_render_hashes() {
  const std::vector<FactoryEntry> &table =
      hs_wasm::get_factory<HASH_W, HASH_H>();
  HS_EXPECT_EQ(table.size(), static_cast<size_t>(HS_EFFECT_COUNT));
  if (!env_var("HS_PRINT_RENDER_HASHES"))
    return;
  for (const FactoryEntry &entry : table) {
    const size_t presets = entry.preset_count > 0 ? entry.preset_count : 1;
    for (size_t preset = 0; preset < presets; ++preset)
      check_preset_render<HASH_W, HASH_H>(entry, preset);
  }
}

/**
 * @brief Module entry point.
 * @return The module's failure count.
 */
inline int run_render_hashes_tests() {
  hs_test::ModuleFixture fixture("render_hashes");
  HS_RUN(test_roster_render_hashes);
  return fixture.result();
}

} // namespace render_hashes_tests
} // namespace hs_test
