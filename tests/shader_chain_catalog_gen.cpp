/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Generator of record for tests/data/shader_chain_catalog.json, the native-ABI
// operator catalog unit_shader_chain golden-pins. Built with the same feature
// set as run_tests so the bytes it writes are the bytes the suite compares.
#include <cstdio>
#include <fstream>
#include <string>

#include "core/engine/engine.h"
#include "core/render/pullback/catalog_export.h"

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: shader_chain_catalog_gen <output-json>\n");
    return 2;
  }

  std::string catalog;
  Pullback::Interp::append_catalog_json(catalog);
  catalog += '\n';

  std::ofstream out(argv[1], std::ios::binary | std::ios::trunc);
  if (!out) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 1;
  }
  out.write(catalog.data(), static_cast<std::streamsize>(catalog.size()));

  // Closed explicitly: the destructor's flush swallows its own failure, and a
  // short write would leave a truncated golden behind an exit status of 0.
  out.close();
  if (!out) {
    std::fprintf(stderr, "cannot write %s\n", argv[1]);
    return 1;
  }
  return 0;
}
