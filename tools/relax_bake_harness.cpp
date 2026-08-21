/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host relax-bake generator. Compiled with HS_RELAX_BAKE_EXTRACT, every
 * SolidBuilder::relax_baked() call reproduces its payload by running exactly
 * `bake.iterations` smoothing steps and logs a RELAX_BAKE block (see
 * core/mesh/solids.h). Running every bake-bearing generator once therefore
 * emits the full asset stream on stdout; tools/relax_bakes.py parses it into
 * core/mesh/relax_bakes_generated.h. Because host relax is deterministic, the
 * emitted bits load unchanged on both host and device.
 *
 * Compiled with HS_RELAX_BAKE_VERIFY instead, the same sweep asserts each
 * re-derivation against the committed payload (the unit_relax_bake_verify gate).
 */
#include <cstdint>
#include <cstdio>
#include "core/mesh/solids.h"

#if defined(HS_RELAX_BAKE_VERIFY)
// relax_baked() steps the registries must reach. The assertions live inside
// relax_baked(), so a sweep that reaches none of them still exits 0; this floor
// is what distinguishes a passing gate from a gate that scored nothing. Raise it
// for a new baked step; lower it only when one is deliberately retired.
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
