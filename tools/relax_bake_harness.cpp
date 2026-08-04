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
 */
#include <cstdint>
#include "core/mesh/solids.h"

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
  return 0;
}
