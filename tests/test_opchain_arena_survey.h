/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// Persistent-arena survey of the whole Islamic registry built op by op. Every
// entry is replayed leg by leg through the same harness the morph suite gates
// with (tests/test_conway_morph.h), including the entries carrying no shipping
// recipe: their chains are transcribed from their generator functions below, so
// the survey reports what each shape would cost if it were built rather than
// generated whole.

#include <algorithm>
#include <array>
#include <cstdio>

#include "tests/test_conway_morph.h"

namespace hs_test {
namespace opchain_arena_survey {

using conway_morph_tests::ChainPeaks;
using conway_morph_tests::ISLAMIC_PERSISTENT_BUDGET;
using conway_morph_tests::ISLAMIC_SCRATCH_A_BUDGET;
using conway_morph_tests::ISLAMIC_SCRATCH_B_BUDGET;
using Solids::Op;
using Solids::OpStep;
using Solids::Recipe;
using Solids::IslamicStarPatterns::D2R;

/** simple_registry seed indices the transcribed chains build from. */
inline constexpr uint8_t SEED_TRUNCATED_OCTAHEDRON = 8;
inline constexpr uint8_t SEED_TRUNCATED_ICOSAHEDRON = 14;
inline constexpr uint8_t SEED_TRUNCATED_ICOSIDODECAHEDRON = 16;
static_assert(std::string_view(
                  Solids::simple_registry[SEED_TRUNCATED_OCTAHEDRON].name) ==
              "truncatedOctahedron");
static_assert(std::string_view(
                  Solids::simple_registry[SEED_TRUNCATED_ICOSAHEDRON].name) ==
              "truncatedIcosahedron");
static_assert(
    std::string_view(
        Solids::simple_registry[SEED_TRUNCATED_ICOSIDODECAHEDRON].name) ==
    "truncatedIcosidodecahedron");

// Chains transcribed from the registry entries that ship no Recipe. SolidBuilder
// defaults are spelled out: relax() is 8 iterations, snub() is t 0.5 / twist 0.
inline constexpr OpStep TICOSA_AMBO_RELAX100_HK54_NEEDLE_STEPS[] = {
    {Op::AMBO},
    {Op::RELAX, 100.0f},
    {Op::HANKIN, 54.0f * D2R},
    {Op::NEEDLE}};
inline constexpr OpStep TICOSA_HK58_CHAMFER63_STEPS[] = {
    {Op::HANKIN, 58.0f * D2R}, {Op::CHAMFER, 0.63f}};
inline constexpr OpStep TICOSA_AMBO_RELAX_TRUNCATE33_HK64_STEPS[] = {
    {Op::AMBO},
    {Op::RELAX, 217.0f},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 64.0f * D2R}};
inline constexpr OpStep DODECA_BEVEL2_RELAX_GYRO_STEPS[] = {
    {Op::BEVEL, 0.2f}, {Op::RELAX, 100.0f}, {Op::GYRO}};
inline constexpr OpStep TICOSIDODECA_BEVEL5_RELAX_HK77_STEPS[] = {
    {Op::BEVEL, 0.5f}, {Op::RELAX, 100.0f}, {Op::HANKIN, 77.0f * D2R}};
inline constexpr OpStep TOCTA_GYRO_KIS_HK17_STEPS[] = {
    {Op::GYRO}, {Op::KIS}, {Op::HANKIN, 17.0f * D2R}};
inline constexpr OpStep TICOSA_AMBO_RELAX_TRUNCATE001_HK59_STEPS[] = {
    {Op::AMBO},
    {Op::RELAX, 8.0f},
    {Op::TRUNCATE, 0.01f},
    {Op::HANKIN, 59.0f * D2R}};
inline constexpr OpStep TICOSA_AMBO_RELAX_TRUNCATE001_HK73_STEPS[] = {
    {Op::AMBO},
    {Op::RELAX, 8.0f},
    {Op::TRUNCATE, 0.01f},
    {Op::HANKIN, 73.0f * D2R}};
inline constexpr OpStep DODECA_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS[] = {
    {Op::HANKIN, 35.0f * D2R}, {Op::AMBO},
    {Op::HANKIN, 62.0f * D2R}, {Op::AMBO},
    {Op::RELAX, 100.0f},       {Op::HANKIN, 42.0f * D2R}};
inline constexpr OpStep TICOSIDODECA_TRUNCATE50D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, 50.0f * D2R}, {Op::AMBO}, {Op::DUAL}};
inline constexpr OpStep TICOSA_HK54_AMBO_HK72_STEPS[] = {
    {Op::HANKIN, 54.0f * D2R}, {Op::AMBO}, {Op::HANKIN, 72.0f * D2R}};
inline constexpr OpStep TICOSA_TRUNCATE50D_AMBO_DUAL_STEPS[] = {
    {Op::TRUNCATE, 50.0f * D2R}, {Op::AMBO}, {Op::DUAL}};
inline constexpr OpStep ICOSA_SNUB_RELAX_TRUNCATE033_HK62_STEPS[] = {
    {Op::SNUB, 0.5f, 0.0f},
    {Op::RELAX, 8.0f},
    {Op::TRUNCATE, 0.33f},
    {Op::HANKIN, 62.0f * D2R}};

/** One transcribed chain, keyed by the registry entry it mirrors. */
struct Transcribed {
  const char *name;
  Recipe recipe;
};

#define HS_TRANSCRIBED(entry, seed, steps)                                     \
  { entry, {seed, steps, static_cast<uint8_t>(std::size(steps))} }

inline constexpr Transcribed TRANSCRIBED_CHAINS[] = {
    HS_TRANSCRIBED("truncatedIcosahedron_ambo_relax100_hk54_needle",
                   SEED_TRUNCATED_ICOSAHEDRON,
                   TICOSA_AMBO_RELAX100_HK54_NEEDLE_STEPS),
    HS_TRANSCRIBED("truncatedIcosahedron_hk58_chamfer63",
                   SEED_TRUNCATED_ICOSAHEDRON, TICOSA_HK58_CHAMFER63_STEPS),
    HS_TRANSCRIBED("truncatedIcosahedron_ambo_relax_truncate33_hk64",
                   SEED_TRUNCATED_ICOSAHEDRON,
                   TICOSA_AMBO_RELAX_TRUNCATE33_HK64_STEPS),
    HS_TRANSCRIBED("dodecahedron_bevel2_relax_gyro", Solids::SEED_DODECAHEDRON,
                   DODECA_BEVEL2_RELAX_GYRO_STEPS),
    HS_TRANSCRIBED("truncatedIcosidodecahedron_bevel5_relax_hk77",
                   SEED_TRUNCATED_ICOSIDODECAHEDRON,
                   TICOSIDODECA_BEVEL5_RELAX_HK77_STEPS),
    HS_TRANSCRIBED("truncatedOctahedron_gyro_kis_hk17",
                   SEED_TRUNCATED_OCTAHEDRON, TOCTA_GYRO_KIS_HK17_STEPS),
    HS_TRANSCRIBED("truncatedIcosahedron_ambo_relax_truncate001_hankin59",
                   SEED_TRUNCATED_ICOSAHEDRON,
                   TICOSA_AMBO_RELAX_TRUNCATE001_HK59_STEPS),
    HS_TRANSCRIBED("truncatedIcosahedron_ambo_relax_truncate001_hankin73",
                   SEED_TRUNCATED_ICOSAHEDRON,
                   TICOSA_AMBO_RELAX_TRUNCATE001_HK73_STEPS),
    HS_TRANSCRIBED("dodecahedron_hk35_ambo_hk62_ambo_relax_hk42",
                   Solids::SEED_DODECAHEDRON,
                   DODECA_HK35_AMBO_HK62_AMBO_RELAX_HK42_STEPS),
    HS_TRANSCRIBED("truncatedIcosidodecahedron_truncate50d_ambo_dual",
                   SEED_TRUNCATED_ICOSIDODECAHEDRON,
                   TICOSIDODECA_TRUNCATE50D_AMBO_DUAL_STEPS),
    HS_TRANSCRIBED("truncatedIcosahedron_hk54_ambo_hk72",
                   SEED_TRUNCATED_ICOSAHEDRON, TICOSA_HK54_AMBO_HK72_STEPS),
    HS_TRANSCRIBED("truncatedIcosahedron_truncate50d_ambo_dual",
                   SEED_TRUNCATED_ICOSAHEDRON,
                   TICOSA_TRUNCATE50D_AMBO_DUAL_STEPS),
    HS_TRANSCRIBED("icosahedron_snub_relax_truncate033_hankin62",
                   Solids::SEED_ICOSAHEDRON,
                   ICOSA_SNUB_RELAX_TRUNCATE033_HK62_STEPS)};

#undef HS_TRANSCRIBED

/**
 * @brief Chain of a registry entry: its shipping recipe, or the transcription
 *        of its generator when it ships none.
 * @param entry Registry entry.
 * @return Pointer to the chain, or nullptr when neither exists.
 */
inline const Recipe *chain_of(const Solids::Entry &entry) {
  if (entry.recipe)
    return entry.recipe;
  for (const Transcribed &t : TRANSCRIBED_CHAINS)
    if (std::string_view(t.name) == entry.name)
      return &t.recipe;
  return nullptr;
}

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
    const Recipe *chain = chain_of(entry);
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
    std::printf("  [survey] %-56s F=%-5zu legs=%zu persistent=%7zu %s%-7zu  "
                "scratch a=%7zu b=%7zu%s\n",
                rows[i].name, p.faces, p.legs, p.persistent,
                fits ? "under " : "OVER +",
                fits ? ISLAMIC_PERSISTENT_BUDGET - p.persistent
                     : p.persistent - ISLAMIC_PERSISTENT_BUDGET,
                p.scratch_a, p.scratch_b,
                p.supported ? "" : "  (unsweepable step)");
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

} // namespace opchain_arena_survey
} // namespace hs_test
