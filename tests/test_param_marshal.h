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
  Timeline<kW>().clear();
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

inline int run_param_marshal_tests() {
  auto scope = hs_test::begin_module("param_marshal");
#define HS_PARAM_ONE(name) check_one<name>(#name);
  HS_EFFECT_LIST(HS_PARAM_ONE)
#undef HS_PARAM_ONE
  return hs_test::end_module(scope);
}

} // namespace param_marshal_tests
} // namespace hs_test
