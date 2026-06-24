/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM parameter-marshaling layer
 * (targets/wasm/param_marshal.h). The bridge exposes two parallel streams to
 * JS — parameter definitions and per-frame values — and the GUI binds them
 * positionally, so a single order mismatch mis-binds every slider. These tests
 * run the marshaling against every registered effect and
 * assert the value stream stays index-aligned with the definition stream, the
 * bool/float distinction is preserved, and a write-by-name round-trips to the
 * same index. The emscripten::val translation in wasm.cpp is a thin shell over
 * this layer and is exercised by the WASM build, not here.
 */
#pragma once

#include "core/effects.h" // HS_EFFECT_LIST roster
#include "core/canvas.h"
#include "core/memory.h"
#include "targets/wasm/param_marshal.h"
#include "tests/test_harness.h"

#include <cstdio>
#include <string_view>
#include <vector>

namespace hs_test {
namespace param_marshal_tests {

constexpr int kW = 288;
constexpr int kH = 144;

/**
 * @brief Marshals one effect through the WASM bridge and asserts the definition
 *        and value streams stay consistent with the source params.
 * @tparam E Effect template, instantiated at the test canvas size kW x kH.
 * @param unnamed Effect name (unused; kept for call-site symmetry with the
 *        HS_EFFECT_LIST macro expansion).
 * @details Checks equal length, index-aligned name/value/type, and a
 *          write-by-name that round-trips to the same index without disturbing
 *          order. This is the core correctness check per effect.
 */
template <template <int, int> class E>
inline bool check_one(const char *) {
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;

  E<kW, kH> effect;
  effect.init();

  std::vector<hs_wasm::ParamView> views;
  std::vector<float> values;
  hs_wasm::collect_param_views(effect, views);
  hs_wasm::fill_param_values(effect, values);

  // Count the source params independently of the marshaled streams.
  size_t n = 0;
  for (const auto &def : effect.getParameters()) {
    (void)def;
    ++n;
  }
  HS_EXPECT_EQ(views.size(), n);
  HS_EXPECT_EQ(values.size(), n);

  // Both streams come from one ordered pass: for every i, name/value/type must
  // match the source param.
  size_t i = 0;
  for (const auto &def : effect.getParameters()) {
    HS_EXPECT(i < views.size(), "marshaled index within range");
    HS_EXPECT_EQ(std::string_view(views[i].name), std::string_view(def.name));
    HS_EXPECT_EQ(views[i].value, values[i]);
    HS_EXPECT_EQ(views[i].is_bool, def.is_bool());
    HS_EXPECT_EQ(views[i].value, def.get());
    ++i;
  }

  // Write an editable float param BY NAME and confirm it reappears at the same
  // index with the order untouched — guards against setParameter landing on the
  // wrong slider. No editable float param -> return false so the caller tallies
  // the skip and roster drift toward such effects stays visible.
  int target = -1;
  for (size_t k = 0; k < views.size(); ++k) {
    const auto &v = views[k];
    if (!v.is_bool && !v.readonly && !v.animated && v.max > v.min) {
      target = static_cast<int>(k);
      break;
    }
  }
  if (target < 0)
    return false;

  const float lo = views[target].min, hi = views[target].max;
  float newv = lo + 0.5f * (hi - lo);
  if (newv == views[target].value)
    newv = lo + 0.25f * (hi - lo);

  HS_EXPECT(effect.updateParameter(views[target].name, newv),
            "updateParameter by name succeeds for an editable param");

  std::vector<hs_wasm::ParamView> views2;
  hs_wasm::collect_param_views(effect, views2);
  HS_EXPECT_EQ(views2.size(), views.size());
  HS_EXPECT_NEAR(views2[target].value, newv, 1e-3f);
  for (size_t k = 0; k < views2.size(); ++k)
    HS_EXPECT_EQ(std::string_view(views2[k].name), std::string_view(views[k].name));
  return true;
}

/**
 * @brief Tallies an effect's parameter count, tracking the roster maximum.
 * @tparam E Effect template, instantiated at the test canvas size kW x kH.
 * @param max_count In/out running maximum; updated if this effect has more
 *        parameters, so the stability pass can size its reserve to the worst
 *        case.
 */
template <template <int, int> class E>
inline void count_one(size_t &max_count) {
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;

  E<kW, kH> effect;
  effect.init();

  size_t n = 0;
  for (const auto &def : effect.getParameters()) {
    (void)def;
    ++n;
  }
  if (n > max_count)
    max_count = n;
}

/**
 * @brief Verifies the per-frame memory-view stability contract: refilling the
 *        reserved stream vectors never reallocates their backing storage.
 * @tparam E Effect template, instantiated at the test canvas size kW x kH.
 * @param views Reusable definition-stream vector, pre-reserved by the caller.
 * @param values Reusable value-stream vector, pre-reserved by the caller.
 * @param view_data Expected backing pointer of @p views (its .data() before
 *        refill).
 * @param view_cap Expected capacity of @p views.
 * @param value_data Expected backing pointer of @p values (its .data() before
 *        refill).
 * @param value_cap Expected capacity of @p values.
 * @details wasm.cpp exposes the value stream to JS as a raw pointer into WASM
 *          linear memory (`paramValues.data()`), read every frame, reusing a
 *          single vector across frames and effect switches. collect_param_views
 *          / fill_param_values promise that — given a caller that reserves
 *          capacity once up front — clear()+push_back never reallocates, so the
 *          exported address stays valid. Dropping the reserve or pushing past it
 *          would silently hand JS a dangling pointer with no crash on the C++
 *          side. This refills the SAME reserved vectors and asserts backing
 *          storage (.data()) and capacity never move.
 */
template <template <int, int> class E>
inline void check_stability_one(std::vector<hs_wasm::ParamView> &views,
                                std::vector<float> &values,
                                const hs_wasm::ParamView *view_data,
                                size_t view_cap, const float *value_data,
                                size_t value_cap) {
  configure_arenas_default();
  Timeline().clear();
  global_timeline_t = 0;

  E<kW, kH> effect;
  effect.init();

  hs_wasm::collect_param_views(effect, views);
  hs_wasm::fill_param_values(effect, values);

  HS_EXPECT(views.data() == view_data,
            "collect_param_views reused the reserved buffer (no realloc)");
  HS_EXPECT(values.data() == value_data,
            "fill_param_values reused the reserved buffer (no realloc)");
  HS_EXPECT_EQ(views.capacity(), view_cap);
  HS_EXPECT_EQ(values.capacity(), value_cap);
}

/**
 * @brief Module entry point: runs the per-effect stream-consistency check
 *        across the whole roster, then the cross-effect memory-stability check.
 * @return The module's failure count.
 */
/**
 * @brief Freezes the effect roster ORDER, not just its count.
 * @details HS_EFFECT_LIST is the single source of truth for the effect ordinal
 *   the WASM factory enumerates and the JS app surfaces (effect-list order, plus
 *   any index-keyed consumer). The startup check (wasm.cpp) and the per-effect
 *   marshaling below both only guarantee the COUNT and within-effect index
 *   alignment; neither notices a reorder. This independent golden list turns any
 *   reorder/insertion/removal into a deliberate, reviewable diff — if it fires,
 *   update kGoldenRoster on purpose to match the new HS_EFFECT_LIST order.
 *   (Sliders bind by parameter name, so a reorder does not mis-bind a slider; it
 *   shifts the effect ordinal, which is what this pins.)
 */
inline void check_roster_order_pinned() {
  // Independent hand-maintained copy of the intended roster order. Must NOT be
  // generated from HS_EFFECT_LIST, or the comparison becomes a tautology.
  static const char *const kGoldenRoster[] = {
      "BZReactionDiffusion", "ChaoticStrings",     "Comets",
      "DistortedRing",       "DreamBalls",         "Dynamo",
      "FlowField",           "Flyby",              "GnomonicStars",
      "GSReactionDiffusion", "HankinSolids",       "HopfFibration",
      "IslamicStars",        "Liquid2D",           "MeshFeedback",
      "MindSplatter",        "MobiusGrid",         "Moire",
      "PetalFlow",           "Raymarch",           "RingShower",
      "RingSpin",            "ShapeShifter",       "SphericalHarmonics",
      "SplineFlow",          "Thrusters",          "Voronoi"};
  // Actual roster, expanded straight from the X-macro source of truth.
  static const char *const kActualRoster[] = {
#define HS_EFFECT_NAME(name) #name,
      HS_EFFECT_LIST(HS_EFFECT_NAME)
#undef HS_EFFECT_NAME
  };
  constexpr size_t kGoldenN = sizeof(kGoldenRoster) / sizeof(kGoldenRoster[0]);
  constexpr size_t kActualN = sizeof(kActualRoster) / sizeof(kActualRoster[0]);
  HS_EXPECT_EQ(kActualN, kGoldenN);
  HS_EXPECT_EQ(static_cast<int>(kActualN), HS_EFFECT_COUNT);
  const size_t n = kActualN < kGoldenN ? kActualN : kGoldenN;
  for (size_t i = 0; i < n; ++i)
    HS_EXPECT_TRUE(std::string_view(kActualRoster[i]) ==
                   std::string_view(kGoldenRoster[i]));
}

inline int run_param_marshal_tests() {
  auto scope = hs_test::begin_module("param_marshal");
  check_roster_order_pinned();
  // Tally how many effects exercised the by-name round-trip; it is skipped for
  // effects with no editable float param. Surface the split and fail if zero.
  int rt_covered = 0, rt_total = 0;
#define HS_PARAM_ONE(name)                                                     \
  do {                                                                         \
    ++rt_total;                                                                \
    if (check_one<name>(#name))                                               \
      ++rt_covered;                                                            \
  } while (0);
  HS_EFFECT_LIST(HS_PARAM_ONE)
#undef HS_PARAM_ONE
  std::printf("  param-marshal by-name round-trip exercised on %d/%d effects "
              "(%d skipped: no editable float param)\n",
              rt_covered, rt_total, rt_total - rt_covered);
  HS_EXPECT(rt_covered > 0,
            "by-name round-trip must run on at least one effect — the roster "
            "drifted to all-non-editable params and the check covers nothing");

  // Size one pair of vectors to the roster's largest parameter set, then marshal
  // every effect through them (the effect-switch path) and confirm the backing
  // storage never reallocates — the memory-view stability the WASM bridge needs.
  size_t max_count = 0;
#define HS_PARAM_COUNT(name) count_one<name>(max_count);
  HS_EFFECT_LIST(HS_PARAM_COUNT)
#undef HS_PARAM_COUNT

  std::vector<hs_wasm::ParamView> views;
  std::vector<float> values;
  views.reserve(max_count);
  values.reserve(max_count);
  const hs_wasm::ParamView *view_data = views.data();
  const float *value_data = values.data();
  const size_t view_cap = views.capacity();
  const size_t value_cap = values.capacity();

#define HS_PARAM_STAB(name)                                                    \
  check_stability_one<name>(views, values, view_data, view_cap, value_data,    \
                            value_cap);
  HS_EFFECT_LIST(HS_PARAM_STAB)
#undef HS_PARAM_STAB

  return hs_test::end_module(scope);
}

} // namespace param_marshal_tests
} // namespace hs_test
