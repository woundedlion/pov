/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file effect_factory.h
 * @brief Per-(W,H) effect factory and resolution dispatch behind the WASM
 *        engine, free of emscripten types.
 *
 * Split out of engine_bindings.h so the native suite can compile and drive it
 * under live HS_CHECK, ASan/UBSan and -O0, the way param_marshal.h and
 * wasm_predicates.h already are.
 */
#pragma once

#include "core/engine/constants.h"
#include "core/control/registry.h"
#include "core/engine/effects.h" // Includes all effect headers (triggers REGISTER_EFFECT)
#include <iterator>              // std::size — X-macro roster tables
#include <string_view>
#include <vector>

/**
 * @brief Builds a concrete factory table from the self-registering entries.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @return Reference to the lazily-built, static per-(W,H) factory table.
 */
template <int W, int H> const std::vector<FactoryEntry> &get_factory() {
  static std::vector<FactoryEntry> table = []() {
    const auto &regs = EffectRegistry::entries();
    std::vector<FactoryEntry> t(regs.size());
    for (size_t i = 0; i < regs.size(); ++i)
      get_fill_fn<W, H>(regs[i])(t[i]);
    return t;
  }();
  return table;
}

/**
 * @brief Looks up an effect by name in the (W,H) factory.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param name Effect name to look up.
 * @return Pointer to the matching entry, or null if the name is unknown (a
 *         typo'd/stale UI string). The entry belongs to the static per-(W,H)
 *         table and stays valid for the module's lifetime.
 * @details Cheap linear scan used by setEffect() to validate a stale/typo'd UI
 *          string BEFORE it tears down the running effect, so an unknown name is
 *          a transactional no-op rather than a blanked engine; the returned
 *          entry then creates the effect without a second lookup.
 */
template <int W, int H>
const FactoryEntry *find_factory_entry(std::string_view name) {
  for (const auto &entry : get_factory<W, H>())
    if (name == entry.name || name == entry.stable_id)
      return &entry;
  return nullptr;
}

// ---------------------------------------------------------------------------
// The (W,H) resolutions the WASM factory can build are the registry's
// HS_RESOLUTIONS rows (core/control/registry.h), so the runtime dispatch
// here and the per-resolution fill functions the registry generates share one
// list and cannot drift: a resolution the registry can build is dispatchable
// here with no second edit. setResolution()/setEffect()/getEffectSizes() all
// expand it via the X-macro below. A new resolution is one edit, in
// HS_RESOLUTIONS (the effect templates must also be instantiable at that
// <W,H>).
// ---------------------------------------------------------------------------

// Pin every resolution row to the MAX_W×MAX_H pixel-buffer bound.
#define X(W, H)                                                                \
  static_assert((W) <= MAX_W && (H) <= MAX_H,                                  \
                "HS_RESOLUTIONS row exceeds the MAX_W×MAX_H pixel buffer");
HS_RESOLUTIONS(X)
#undef X

/** @brief One HS_RESOLUTIONS row as runtime values. */
struct WasmResolution {
  int w; /**< Canvas width in pixels. */
  int h; /**< Canvas height in pixels. */
};

// The rows as data, so the constructor can bootstrap on the first one instead of
// naming a preset that HS_RESOLUTIONS may later drop.
inline constexpr WasmResolution WASM_RESOLUTIONS[] = {
#define X(W, H) {(W), (H)},
    HS_RESOLUTIONS(X)
#undef X
};
static_assert(std::size(WASM_RESOLUTIONS) > 0,
              "HS_RESOLUTIONS must carry a row to bootstrap on");

// The registered effect names, in HS_EFFECT_LIST order. Only the first is used
// (the constructor's bootstrap), for the same reason as WASM_RESOLUTIONS.
inline constexpr const char *WASM_EFFECT_NAMES[] = {
#define X(name) #name,
    HS_EFFECT_LIST(X)
#undef X
};
static_assert(std::size(WASM_EFFECT_NAMES) == HS_EFFECT_COUNT,
              "HS_EFFECT_LIST and HS_EFFECT_COUNT must agree");

/**
 * @brief Invokes f.operator()<W,H>() for the single HS_RESOLUTIONS row
 *        matching the runtime (w,h).
 * @tparam F Type of the templated callable.
 * @param w Runtime canvas width in pixels to match against the row list.
 * @param h Runtime canvas height in pixels to match against the row list.
 * @param f A C++20 templated callable (e.g. `[]<int W, int H>(){...}`) run with
 *          the matching row's dimensions as compile-time template arguments.
 * @return true iff a row matched and f was invoked; false otherwise.
 * @details One shared expansion of the resolution list for every runtime
 *          dispatch site (setEffect validate/create, resolution checks), so they
 *          can never drift and a new per-resolution step lands in exactly one
 *          place.
 */
template <typename F> inline bool dispatch_resolution(int w, int h, F &&f) {
#define X(W, H)                                                                \
  if (w == (W) && h == (H)) {                                                  \
    f.template operator()<(W), (H)>();                                         \
    return true;                                                               \
  }
  HS_RESOLUTIONS(X)
#undef X
  return false;
}

/**
 * @brief Reports whether (w,h) is a resolution the factory can build.
 * @param w Candidate canvas width in pixels.
 * @param h Candidate canvas height in pixels.
 * @return true iff (w,h) is one of the HS_RESOLUTIONS rows.
 */
inline bool wasm_resolution_supported(int w, int h) {
  return dispatch_resolution(w, h, []<int W, int H>() {});
}
