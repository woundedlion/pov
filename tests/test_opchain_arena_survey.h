/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

// Persistent-arena survey of the whole Islamic registry built op by op. Every
// entry is replayed leg by leg through the same harness the morph suite gates
// with (tests/test_conway_morph.h), off the shipping Recipe that
// test_solids.h's bitwise gate pins to the entry's generator.

#include <algorithm>
#include <array>
#include <cstdio>

#include "tests/test_conway_morph.h"

namespace hs_test {
namespace opchain_arena_survey_tests {

using conway_morph_tests::ChainPeaks;
using conway_morph_tests::ISLAMIC_PERSISTENT_BUDGET;
using conway_morph_tests::ISLAMIC_SCRATCH_A_BUDGET;
using conway_morph_tests::ISLAMIC_SCRATCH_B_BUDGET;
using Solids::Recipe;

/** One surveyed row. */
struct Row {
  const char *name;
  ChainPeaks peaks;
};

/**
 * @brief Replays every Islamic registry entry as a build chain and reports the
 *        persistent and scratch high-waters against IslamicStars' split.
 */
inline void test_islamic_registry_arena_survey() {
  const std::span<const Solids::Entry> entries =
      Solids::Collections::get_islamic_solids();
  std::array<Row, std::size(Solids::islamic_registry)> rows{};
  size_t n = 0;

  for (const Solids::Entry &entry : entries) {
    const Recipe *chain = entry.recipe;
    HS_EXPECT_TRUE(chain != nullptr);
    if (!chain)
      continue;
    rows[n].name = entry.name;
    rows[n].peaks = conway_morph_tests::replay_build_chain(entry.name, *chain,
                                                           /*gate=*/false);
    ++n;
  }
  HS_EXPECT_EQ(n, entries.size());

  std::sort(rows.begin(), rows.begin() + n, [](const Row &x, const Row &y) {
    return x.peaks.persistent > y.peaks.persistent;
  });

  size_t over = 0, worst_a = 0, worst_b = 0;
  std::printf("  [survey] persistent budget %zu B (scratch a %zu B, b %zu B)\n",
              (size_t)ISLAMIC_PERSISTENT_BUDGET,
              (size_t)ISLAMIC_SCRATCH_A_BUDGET,
              (size_t)ISLAMIC_SCRATCH_B_BUDGET);
  for (size_t i = 0; i < n; ++i) {
    const ChainPeaks &p = rows[i].peaks;
    const bool fits = p.persistent <= ISLAMIC_PERSISTENT_BUDGET;
    over += fits ? 0 : 1;
    worst_a = std::max(worst_a, p.scratch_a);
    worst_b = std::max(worst_b, p.scratch_b);
    std::printf(
        "  [survey] %-56s F=%-5zu legs=%zu persistent=%7zu %s%-7zu  "
        "scratch a=%7zu b=%7zu%s\n",
        rows[i].name, p.faces, p.legs, p.persistent, fits ? "under " : "OVER +",
        fits ? ISLAMIC_PERSISTENT_BUDGET - p.persistent
             : p.persistent - ISLAMIC_PERSISTENT_BUDGET,
        p.scratch_a, p.scratch_b, p.supported ? "" : "  (unsweepable step)");
  }
  std::printf("  [survey] %zu of %zu entries fit; worst scratch a=%zu / %zu, "
              "b=%zu / %zu\n",
              n - over, n, worst_a, (size_t)ISLAMIC_SCRATCH_A_BUDGET, worst_b,
              (size_t)ISLAMIC_SCRATCH_B_BUDGET);
  // Every recipe now fits the persistent budget; gating at zero keeps a
  // build-path regression from silently pushing any shape past it.
  HS_EXPECT_LE(over, (size_t)0);
  HS_EXPECT_LE(worst_a, ISLAMIC_SCRATCH_A_BUDGET);
  HS_EXPECT_LE(worst_b, ISLAMIC_SCRATCH_B_BUDGET);
}

/**
 * @brief Runs the opchain arena survey.
 * @return The module's failure count.
 */
inline int run_opchain_arena_survey_tests() {
  hs_test::ModuleFixture fixture("opchain_arena_survey");
  test_islamic_registry_arena_survey();
  return fixture.result();
}

} // namespace opchain_arena_survey_tests
} // namespace hs_test
