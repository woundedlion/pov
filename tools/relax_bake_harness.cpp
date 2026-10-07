/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host relax-bake generator. With HS_RELAX_BAKE_EXTRACT, each
 * SolidBuilder::relax_baked() call re-derives its payload with
 * `bake.iterations` as the smoothing cap and logs a RELAX_BAKE block on stdout
 * for tools/relax_bakes.py. With HS_RELAX_BAKE_VERIFY, each re-derivation is
 * asserted against the committed payload.
 *
 * New bakes go in core/mesh/relax_bake_specs.h; main() must reach their
 * generator, and MIN_RELAX_BAKES_VERIFIED rises by one per added step.
 */
#include <cstdint>
#include <cstdio>
#include "core/mesh/solids.h"

#if defined(HS_RELAX_BAKE_VERIFY)
// Floor on relax_baked() steps reached; without it a sweep reaching none passes.
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

  // Shared payloads re-emit identically; relax_bakes.py de-duplicates by name.
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
