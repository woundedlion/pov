/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM engine's factory layer
 * (targets/wasm/effect_factory.h): the per-(W,H) FactoryEntry tables built from
 * the self-registering effects, the name lookup setEffect() validates against,
 * and the HS_RESOLUTIONS runtime dispatch. FactoryEntry::creator is the render
 * path's one std::make_unique; the browser is otherwise its only caller, at
 * -O3 -flto with assertions off, so it is driven here under live HS_CHECK and
 * the sanitizers.
 */
#pragma once

#include <cstdio>
#include <cstring>
#include <iterator>
#include <memory>
#include <string_view>
#include "core/render/canvas.h"
#include "targets/wasm/effect_factory.h"
#include "tests/test_effects.h" // reset_effect_globals
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace effect_factory_tests {

/** @brief Frames rendered per factory-built effect. */
constexpr int FACTORY_FRAMES = 2;

/**
 * @brief Checks one resolution's factory table against the registry.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Every registration must yield one entry with a name, a creator and a
 *          non-zero size, and the table must be the same object on every call —
 *          setEffect() hands out pointers into it and holds them across frames.
 */
template <int W, int H> inline void verify_factory_table() {
  const std::vector<FactoryEntry> &table = get_factory<W, H>();
  const std::vector<FactoryEntry> &again = get_factory<W, H>();
  HS_EXPECT_EQ(table.size(), EffectRegistry::entries().size());
  HS_EXPECT_EQ(table.size(), static_cast<size_t>(HS_EFFECT_COUNT));
  HS_EXPECT_TRUE(&table == &again);

  for (const FactoryEntry &entry : table) {
    HS_EXPECT_TRUE(!entry.name.empty());
    HS_EXPECT_TRUE(static_cast<bool>(entry.creator));
    HS_EXPECT_TRUE(entry.size > 0);
  }
}

/**
 * @brief Checks the name lookup setEffect() validates a UI string with.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Every HS_EFFECT_LIST name must resolve to the entry carrying it, and
 *          an unregistered name must return null rather than a neighbouring
 *          entry — setEffect() tears the running effect down only after this
 *          lookup succeeds.
 */
template <int W, int H> inline void verify_factory_lookup() {
  const auto lookup = [](std::string_view name) {
    return find_factory_entry<W, H>(name);
  };
  for (const char *name : WASM_EFFECT_NAMES) {
    const FactoryEntry *entry = lookup(name);
    HS_EXPECT_TRUE(entry != nullptr);
    if (entry)
      HS_EXPECT_TRUE(entry->name == name);
  }
  HS_EXPECT_TRUE(lookup("") == nullptr);
  HS_EXPECT_TRUE(lookup("NoSuchEffect") == nullptr);
  // A prefix of a registered name must miss: the scan compares whole names.
  const std::string_view first = WASM_EFFECT_NAMES[0];
  HS_EXPECT_TRUE(lookup(first.substr(0, first.size() - 1)) == nullptr);
}

/**
 * @brief Builds, renders and destroys every registered effect at one resolution.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The heap-allocated lifetime the browser drives: creator() ->
 *          init() -> frames -> reset(). Each instance is destroyed before the
 *          next is built, since Effect permits one live instance at a time, and
 *          the globals are reset per effect the way the roster sweeps do.
 */
template <int W, int H> inline void drive_factory_lifecycle() {
  std::printf("  -- factory create/render/destroy @ %dx%d --\n", W, H);
  for (const FactoryEntry &entry : get_factory<W, H>()) {
    effects_tests::reset_effect_globals();

    std::unique_ptr<Effect> effect = entry.creator();
    HS_EXPECT_TRUE(effect != nullptr);
    if (!effect)
      continue;
    HS_EXPECT_EQ(effect->width(), W);
    HS_EXPECT_EQ(effect->height(), H);

    effect->init();
    for (int f = 0; f < FACTORY_FRAMES; ++f) {
      effect->draw_frame();
      // Consume the queued frame, else the next Canvas ctor spin-waits.
      effect->advance_display();
    }
    effect.reset();
  }
}

/**
 * @brief Exercises the HS_RESOLUTIONS runtime dispatch.
 * @details Every declared row must dispatch and carry its own (W,H) into the
 *          callable; a size off the list must dispatch nowhere rather than fall
 *          into the nearest row.
 */
inline void test_resolution_dispatch() {
  for (const WasmResolution &row : WASM_RESOLUTIONS) {
    HS_EXPECT_TRUE(wasm_resolution_supported(row.w, row.h));
    int seen_w = 0, seen_h = 0, hits = 0;
    const bool dispatched =
        dispatch_resolution(row.w, row.h, [&]<int W, int H>() {
          seen_w = W;
          seen_h = H;
          ++hits;
        });
    HS_EXPECT_TRUE(dispatched);
    HS_EXPECT_EQ(hits, 1);
    HS_EXPECT_EQ(seen_w, row.w);
    HS_EXPECT_EQ(seen_h, row.h);
  }

  HS_EXPECT_TRUE(!wasm_resolution_supported(0, 0));
  HS_EXPECT_TRUE(!wasm_resolution_supported(-1, -1));
  HS_EXPECT_TRUE(!wasm_resolution_supported(WASM_RESOLUTIONS[0].w,
                                            WASM_RESOLUTIONS[0].h + 1));
  int spurious = 0;
  const bool off_list =
      dispatch_resolution(97, 21, [&]<int W, int H>() { ++spurious; });
  HS_EXPECT_TRUE(!off_list);
  HS_EXPECT_EQ(spurious, 0);
}

/**
 * @brief Checks the roster tables the engine bootstraps from.
 * @details The constructor starts on WASM_RESOLUTIONS[0] / WASM_EFFECT_NAMES[0]
 *          rather than a named preset, so both must be non-empty and the first
 *          row buildable.
 */
inline void test_bootstrap_rows() {
  HS_EXPECT_TRUE(std::size(WASM_RESOLUTIONS) > 0);
  HS_EXPECT_EQ(std::size(WASM_EFFECT_NAMES),
               static_cast<size_t>(HS_EFFECT_COUNT));
  HS_EXPECT_TRUE(
      wasm_resolution_supported(WASM_RESOLUTIONS[0].w, WASM_RESOLUTIONS[0].h));
  HS_EXPECT_TRUE(std::strlen(WASM_EFFECT_NAMES[0]) > 0);
}

/**
 * @brief Checks both resolutions' factory tables and name lookups.
 */
inline void test_factory_tables() {
  verify_factory_table<96, 20>();
  verify_factory_table<288, 144>();
  verify_factory_lookup<96, 20>();
  verify_factory_lookup<288, 144>();
}

/**
 * @brief Drives the create/render/destroy lifecycle at every resolution.
 * @details The 288x144 leg carries the FULL tier's cost (a heap instance of
 *          every effect at the production resolution), so it runs where the
 *          roster's other full-resolution sweeps do.
 */
inline void test_factory_lifecycle() {
  drive_factory_lifecycle<96, 20>();
  if (effects_tests::effects_full_suite())
    drive_factory_lifecycle<288, 144>();
}

/**
 * @brief Module entry point: runs the factory and dispatch checks.
 * @return The module's failure count.
 */
inline int run_effect_factory_tests() {
  hs_test::ModuleFixture fixture("effect_factory");
  test_bootstrap_rows();
  test_resolution_dispatch();
  test_factory_tables();
  test_factory_lifecycle();
  return fixture.result();
}

} // namespace effect_factory_tests
} // namespace hs_test
