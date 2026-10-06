/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host relax-bake generator. Compiled with HS_RELAX_BAKE_EXTRACT, every
 * SolidBuilder::relax_baked() call reproduces its payload by using
 * `bake.iterations` as the smoothing iteration cap and logs a RELAX_BAKE block,
 * so running every bake-bearing generator once emits the full asset stream on
 * stdout for tools/relax_bakes.py.
 *
 * Authoring: add or retune names and iterations in core/mesh/relax_bake_specs.h.
 * Ensure main() reaches each new generator, rebuild relax_bake_gen, regenerate,
 * and run relax_bake_verify. Raise MIN_RELAX_BAKES_VERIFIED for each added step.
 *
 * Compiled with HS_RELAX_BAKE_VERIFY instead, the same sweep asserts each
 * re-derivation against the committed payload.
 */
#include <cstdint>
#include <cstdio>
#include "core/mesh/solids.h"

#if defined(HS_RELAX_BAKE_VERIFY)
// relax_baked() steps the registries must reach; a sweep reaching none of them
// still exits 0.
static constexpr int MIN_RELAX_BAKES_VERIFIED = 21;
#endif

int main() {
  static uint8_t arena_a[1 << 22];
  static uint8_t arena_b[1 << 22];

  auto run = [&](const Solids::Entry &e) {
    Arena a(arena_a, sizeof(arena_a));
    Arena b(arena_b, sizeof(arena_b));
    e.generate(a, b);
  };

  // Duplicate payloads (the ambo prefix is shared by several stars) re-emit
  // identically and are de-duplicated by name downstream.
  for (auto reg : Solids::all_registries())
    for (const auto &e : reg)
      run(e);
#if defined(HS_RELAX_BAKE_VERIFY)
  if (Solids::relax_bakes_verified < MIN_RELAX_BAKES_VERIFIED) {
    std::printf("relax bake verify: reached %d relax_baked() steps of %d — a "
                "recipe dropped one, so its committed payload is pinned to "
                "nothing; see tools/relax_bake_harness.cpp\n",
                Solids::relax_bakes_verified, MIN_RELAX_BAKES_VERIFIED);
    return 1;
  }
  std::printf("relax bake verify: %d relax_baked() steps re-derived\n",
              Solids::relax_bakes_verified);
#endif
  return 0;
}
