/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Roster-wide effect sweeps: smoke, cross-run determinism, clip-clear parity
 * and paused rendering. Assertions are same-binary comparisons or coarse
 * properties with no fixed IEEE reference values, so the module runs on the
 * -ffast-math -fno-finite-math-only axis.
 */
#pragma once

#include "tests/test_effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include <cstdlib>

namespace hs_test {
namespace effects_smoke_tests {

using namespace hs_test::effects_tests;

constexpr int PAUSED_FRAMES = 4;

/** @brief The roster slice one process sweeps at the production resolution. */
struct EffectShard {
  int index = 0;        ///< This shard's index, in [0, count).
  int count = 1;        ///< Number of shards the roster is split across.
  bool sharded = false; ///< Whether HS_EFFECTS_SHARD selected a slice.
  bool valid = true;    ///< False when HS_EFFECTS_SHARD is malformed.

  /**
   * @brief Whether this shard sweeps a roster entry.
   * @param roster_index Entry's position in HS_EFFECT_LIST.
   * @return True for entries whose index % count == index.
   */
  bool takes(int roster_index) const { return roster_index % count == index; }
};

/**
 * @brief Reads the roster slice from HS_EFFECTS_SHARD.
 * @return Every effect when unset or empty; "k/N" selects shard k of N.
 */
inline EffectShard effect_shard_from_env() {
  EffectShard shard;
  const char *spec = hs_test::env_var("HS_EFFECTS_SHARD");
  if (!spec || spec[0] == '\0')
    return shard;
  shard.sharded = true;
  char *end = nullptr;
  const long index = std::strtol(spec, &end, 10);
  if (end == spec || *end != '/') {
    shard.valid = false;
    return shard;
  }
  const char *count_text = end + 1;
  const long count = std::strtol(count_text, &end, 10);
  shard.valid = end != count_text && *end == '\0' && count >= 1 &&
                count <= 64 && index >= 0 && index < count;
  if (shard.valid) {
    shard.index = static_cast<int>(index);
    shard.count = static_cast<int>(count);
  }
  return shard;
}

template <template <int, int> class E, int W = SMALL_W, int H = SMALL_H>
inline void paused_render_one(const char *name) {
  HS_CONTEXT(name);
  reset_globals();
  pin_frame_clock(0);

  const int FRAMES = PAUSED_FRAMES;
  E<W, H> effect;
  effect.setAnimationsPaused(true);
  effect.init();
  HS_EXPECT_TRUE(effect.animations_paused());
  for (int f = 0; f < FRAMES; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
  }

  const uint64_t acc = frame_energy<W, H>(effect);

  if (acc == 0)
    std::printf(
        "  PAUSED-BLANK %-20s final paused frame (of %d) has no lit pixel "
        "@ %dx%d\n",
        name, FRAMES, W, H);
  HS_EXPECT(acc > 0, "effect must produce non-black output while paused");
}

/**
 * @brief Each roster effect produces visible output while its presets are
 * paused.
 */
inline void test_every_effect_renders_while_paused() {
  std::printf("  -- paused render, %d frames --\n", PAUSED_FRAMES);
#define HS_PAUSED_ONE(name) paused_render_one<name>(#name);
  HS_EFFECT_LIST(HS_PAUSED_ONE)
#undef HS_PAUSED_ONE
}

/**
 * @brief Module entry point for the roster-wide effect sweeps.
 * @return Module result code from hs_test::end_module (0 on success).
 * @details HS_EFFECTS_FULL=1 adds the smoke, determinism and paused-render
 * passes at the production resolution. HS_EFFECTS_SHARD=k/N runs only those
 * passes, for roster shard k of N.
 */
inline int run_effects_smoke_tests() {
  hs_test::ModuleFixture fixture("effects_smoke");
  const EffectShard shard = effect_shard_from_env();
  HS_EXPECT_TRUE(shard.valid);
  HS_EXPECT_TRUE(!shard.sharded || effects_full_suite());
  int roster_index = 0;

  if (effects_full_suite()) {
    // Full production-resolution roster passes: smoke, then cross-run
    // determinism under the injected clock.
#define HS_SMOKE_ONE(name)                                                     \
  if (shard.takes(roster_index++))                                             \
    smoke_one<name>(#name);
    HS_EFFECT_LIST(HS_SMOKE_ONE)
#undef HS_SMOKE_ONE
    roster_index = 0;
#define HS_DET_ONE(name)                                                       \
  if (shard.takes(roster_index++))                                             \
    determinism_one<name>(#name);
    HS_EFFECT_LIST(HS_DET_ONE)
#undef HS_DET_ONE
  } else {
    std::printf("  [TIER] production-resolution smoke and determinism omitted; "
                "set HS_EFFECTS_FULL=1\n");
  }
  if (shard.sharded) {
    roster_index = 0;
#define HS_PAUSED_SHARD(name)                                                  \
  if (shard.takes(roster_index++))                                             \
    paused_render_one<name, DEFAULT_W, DEFAULT_H>(#name);
    HS_EFFECT_LIST(HS_PAUSED_SHARD)
#undef HS_PAUSED_SHARD
    return fixture.result();
  }

  // Small-aspect <96,20> roster passes run in every tier.
  std::printf("  -- small-aspect resolution %dx%d --\n", SMALL_W, SMALL_H);
#define HS_SMOKE_ONE_SMALL(name) smoke_one<name, SMALL_W, SMALL_H>(#name);
  HS_EFFECT_LIST(HS_SMOKE_ONE_SMALL)
#undef HS_SMOKE_ONE_SMALL
#define HS_DET_ONE_SMALL(name) determinism_one<name, SMALL_W, SMALL_H>(#name);
  HS_EFFECT_LIST(HS_DET_ONE_SMALL)
#undef HS_DET_ONE_SMALL

  // Check clip-clear parity across the effect roster.
  std::printf("  -- clip-clear display parity --\n");
#define HS_CLIP_PARITY_ONE(name)                                               \
  clip_clear_parity_one<name, SMALL_W, SMALL_H>(#name);
  HS_EFFECT_LIST(HS_CLIP_PARITY_ONE)
#undef HS_CLIP_PARITY_ONE

  HS_RUN(test_every_effect_renders_while_paused);
  if (effects_full_suite()) {
#define HS_PAUSED_FULL(name)                                                   \
  paused_render_one<name, DEFAULT_W, DEFAULT_H>(#name);
    HS_EFFECT_LIST(HS_PAUSED_FULL)
#undef HS_PAUSED_FULL
  } else {
    std::printf(
        "  [TIER] production-resolution paused renders omitted; set HS_EFFECTS_FULL=1\n");
  }

  return fixture.result();
}

} // namespace effects_smoke_tests
} // namespace hs_test
