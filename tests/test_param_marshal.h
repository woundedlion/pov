/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM parameter-marshaling layer
 * (targets/wasm/param_marshal.h). The bridge exposes two parallel streams to
 * JS — parameter definitions and per-frame values — and the GUI binds them
 * positionally, so a single order mismatch mis-binds every slider (the P2-12
 * class). These tests run the marshaling against every registered effect and
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

#include <vector>

namespace hs_test {
namespace param_marshal_tests {

constexpr int kW = 288;
constexpr int kH = 144;

template <template <int, int> class E>
inline void check_one(const char *name) {
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

  // Definition stream and value stream are derived from one ordered pass, so
  // for every i: same name, same value, same type as the source param.
  size_t i = 0;
  for (const auto &def : effect.getParameters()) {
    HS_EXPECT(i < views.size(), "marshaled index within range");
    HS_EXPECT_EQ(views[i].name, def.name);     // order matches the source
    HS_EXPECT_EQ(views[i].value, values[i]);   // value stream index-aligned
    HS_EXPECT_EQ(views[i].is_bool, def.is_bool());
    HS_EXPECT_EQ(views[i].value, def.get());
    ++i;
  }

  // Write an editable float param BY NAME and confirm it reappears at the same
  // index with the rest of the order untouched — the exact mis-binding P2-12
  // warned about (a value pushed back through setParameter landing on the wrong
  // slider). Skip effects with no editable float param.
  int target = -1;
  for (size_t k = 0; k < views.size(); ++k) {
    const auto &v = views[k];
    if (!v.is_bool && !v.readonly && !v.animated && v.max > v.min) {
      target = static_cast<int>(k);
      break;
    }
  }
  if (target >= 0) {
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
      HS_EXPECT_EQ(views2[k].name, views[k].name); // order/names stable
  }
}

// Tally an effect's parameter count so the stability pass below can size its
// reserve to the roster's worst case.
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

// Memory-stability contract. wasm.cpp exposes the value stream to JS as a raw
// pointer into WASM linear memory (`paramValues.data()`), read every frame, and
// reuses a single vector across frames AND effect switches. collect_param_views
// / fill_param_values promise that — given a caller that reserves capacity once
// up front — clear()+push_back never reallocates, so the exported address stays
// valid. A regression that dropped the reserve or pushed past it would silently
// hand JS a dangling pointer with no crash on the C++ side. Refill the SAME
// reserved vector and assert its backing storage (.data()) and capacity never
// move.
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

inline int run_param_marshal_tests() {
  auto scope = hs_test::begin_module("param_marshal");
#define HS_PARAM_ONE(name) check_one<name>(#name);
  HS_EFFECT_LIST(HS_PARAM_ONE)
#undef HS_PARAM_ONE

  // Size a single pair of vectors to the roster's largest parameter set, then
  // marshal every effect through them in turn (the effect-switch path) and
  // confirm the backing storage is never reallocated — the per-frame memory-view
  // stability the WASM bridge depends on.
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
