/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM factory layer (targets/wasm/effect_factory.h).
 */
#pragma once

#include <cstdio>
#include <limits>
#include <vector>
#include <memory>
#include <string_view>
#include "core/render/canvas.h"
#include "targets/Phantasm/phantasm_playlist.h"
#include "targets/wasm/effect_factory.h"
#include "tests/test_effects.h" // reset_effect_globals
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace effect_factory_tests {

/** @brief Minimum frames rendered per factory-built effect. */
constexpr int FACTORY_FRAMES = 2;
/** @brief Maximum startup window for a factory-built effect to produce output. */
constexpr int FACTORY_OUTPUT_WINDOW = 64;

/** @brief Frozen compatibility keys and seed hashes for installed effects. */
struct IdentityGolden {
  std::string_view name;     ///< Class/browser lookup alias.
  std::string_view seed_key; ///< Immutable RNG compatibility key.
  uint64_t seed;             ///< Frozen 64-bit seed hash.
};
/// Original effect RNG streams, independent of the identity roster.
constexpr std::array IDENTITY_GOLDENS{
    IdentityGolden{"AlienBrain", "alien-brain", 0x1B7D2D407B3D2D0DULL},
    IdentityGolden{"AlienCore", "alien-core", 0xD295C46DF9FAC866ULL},
    IdentityGolden{"AlienOcean", "alien-ocean", 0x210664B0D66AFB48ULL},
    IdentityGolden{"AshCloud", "ash-cloud", 0xC252C918EF8F2C62ULL},
    IdentityGolden{"BZReactionDiffusion", "BZReactionDiffusion",
                   0xEA36F0CB4D53FE03ULL},
    IdentityGolden{"ChromaticLichen", "chromatic-lichen",
                   0xBDB34E556EC2CB67ULL},
    IdentityGolden{"Comets", "Comets", 0x20CC05B1BB853C06ULL},
    IdentityGolden{"CosmicEyeball", "cosmic-eyeball", 0xC66D1C858490A981ULL},
    IdentityGolden{"DisplacementField", "DisplacementField",
                   0x66DEEE04DCD50C4FULL},
    IdentityGolden{"DreamBalls", "DreamBalls", 0xB9C2D5D2EF43B3DAULL},
    IdentityGolden{"Dynamo", "Dynamo", 0xE278943A66B53DF4ULL},
    IdentityGolden{"Fishbowl", "Fishbowl", 0x3F899529C749D7A0ULL},
    IdentityGolden{"Galaxies", "Galaxies", 0x9874FB9E47360D77ULL},
    IdentityGolden{"GnomonicStars", "GnomonicStars", 0x025960528523F2FEULL},
    IdentityGolden{"GridSpace", "grid-space", 0x816DB39119548BD2ULL},
    IdentityGolden{"GSReactionDiffusion", "GSReactionDiffusion",
                   0xFA20DF830B51B72FULL},
    IdentityGolden{"HankinSolids", "HankinSolids", 0x30B3BDF0CC9BF813ULL},
    IdentityGolden{"HopfFibration", "HopfFibration", 0x4D6A1E9746B48E48ULL},
    IdentityGolden{"HyperLattice", "HyperLattice", 0x6954B5B68BA922CEULL},
    IdentityGolden{"IslamicStars", "IslamicStars", 0xB085099D2EA92F09ULL},
    IdentityGolden{"JewelMelt", "jewel-melt", 0x59D1DF6060825FE6ULL},
    IdentityGolden{"KaleidoscopeFlowers", "kaleidoscope-flowers",
                   0x00B1D40F7728DCDCULL},
    IdentityGolden{"KaleidoscopeHexBright", "kaleidoscope-hex-bright",
                   0xFD71C050E167C355ULL},
    IdentityGolden{"KaleidoscopeHexOil", "kaleidoscope-hex-oil",
                   0x0DB2472904AF2C50ULL},
    IdentityGolden{"KaleidoscopeHexSoft", "kaleidoscope-hex-soft",
                   0x2FA8CC12BC456AEAULL},
    IdentityGolden{"KaleidoscopeMandala", "kaleidoscope-mandala",
                   0xF65DEA98AF8BFB80ULL},
    IdentityGolden{"KaleidoscopePentBright", "kaleidoscope-pent-bright",
                   0xAE0AE13EF286A15BULL},
    IdentityGolden{"KaleidoscopeSmooth", "kaleidoscope-smooth",
                   0xDC60947F5EF28B5CULL},
    IdentityGolden{"KaleidoscopeStainedGlass", "kaleidoscope-stained-glass",
                   0x458439924D917BF4ULL},
    IdentityGolden{"LatticeMelt", "lattice-melt", 0x8A8239B618698E89ULL},
    IdentityGolden{"MermaidSkin", "mermaid-skin", 0xF49045D62E42DE36ULL},
    IdentityGolden{"MeshFeedback", "MeshFeedback", 0xC445C3A1098582B6ULL},
    IdentityGolden{"MindSplatter", "MindSplatter", 0x206A75CB20AE4025ULL},
    IdentityGolden{"MobiusGrid", "mobius-grid", 0x090E86805292EE94ULL},
    IdentityGolden{"MobiusRings", "MobiusRings", 0xB7C98BEE3FB3F04EULL},
    IdentityGolden{"PetalFlow", "PetalFlow", 0xC3F812CFF88DFE2AULL},
    IdentityGolden{"Raymarch", "Raymarch", 0x64165ACC0D1BEA95ULL},
    IdentityGolden{"RingShower", "RingShower", 0x8A0968F09B972CE9ULL},
    IdentityGolden{"RingSpin", "RingSpin", 0x128AE1065F9277BBULL},
    IdentityGolden{"ShapeShifter", "ShapeShifter", 0xE50A50290CAC42C0ULL},
    IdentityGolden{"SphericalHarmonics", "SphericalHarmonics",
                   0x9CB01486B8230EFEULL},
    IdentityGolden{"Thrusters", "Thrusters", 0x61547828E9116CB4ULL},
    IdentityGolden{"Voronoi", "Voronoi", 0xBC5209F3455BC431ULL},
    IdentityGolden{"ShaderChain", "shader-chain", 0x1968066A5CD68D2BULL},
};

/**
 * @brief Finds a frozen RNG identity by class/browser name.
 * @param name Class/browser name.
 * @return Matching snapshot, or null for an unknown effect.
 */
constexpr const IdentityGolden *identity_golden(std::string_view name) {
  for (const auto &value : IDENTITY_GOLDENS)
    if (value.name == name)
      return &value;
  return nullptr;
}

/** @brief Pins canonical ID types and immutable aliases/seeds for every effect. */
inline void test_identity_contract() {
#define HS_VERIFY_EFFECT_IDENTITY(Name, Id, Display, Legacy, Seed)             \
  do {                                                                         \
    using E = Name<96, 20>;                                                    \
    static_assert(                                                             \
        std::is_same_v<decltype(E::EFFECT_ID), const std::string_view>);       \
    constexpr auto identity = E::EFFECT_IDENTITY;                              \
    HS_EXPECT_TRUE(hs::canonical_effect_id(identity.machine_id));              \
    HS_EXPECT_TRUE(identity.machine_id == E::EFFECT_ID);                       \
    const auto *golden = identity_golden(#Name);                               \
    HS_EXPECT_TRUE(golden != nullptr);                                         \
    if (golden) {                                                              \
      HS_EXPECT_TRUE(identity.display_name == golden->name);                   \
      HS_EXPECT_TRUE(identity.legacy_id == golden->seed_key);                  \
      HS_EXPECT_TRUE(identity.rng_seed_key == golden->seed_key);               \
      HS_EXPECT_EQ(hs::stable_effect_seed(identity.rng_seed_key),              \
                   golden->seed);                                              \
    }                                                                          \
  } while (false);
  HS_EFFECT_IDENTITY_LIST(HS_VERIFY_EFFECT_IDENTITY)
#undef HS_VERIFY_EFFECT_IDENTITY
  constexpr hs::EffectIdentity ALPHA{"alpha", "Alpha", "old-alpha",
                                     "rng-alpha"};
  static_assert(ALPHA.matches("alpha") && ALPHA.matches("Alpha") &&
                ALPHA.matches("old-alpha"));
  static_assert(!ALPHA.matches("") && !ALPHA.matches("ALPHA") &&
                !ALPHA.matches("alph"));
  EffectRegistration first{};
  first.identity = ALPHA;
  EffectRegistration second{};
  second.identity = {"beta", "Beta", "old-beta", "rng-beta"};
  HS_EXPECT_TRUE(registration_names_unique(std::array{first, second}));
  second.identity.legacy_id = "old-alpha";
  HS_EXPECT_FALSE(registration_names_unique(std::array{first, second}));
  second.identity.legacy_id = "alpha";
  HS_EXPECT_FALSE(registration_names_unique(std::array{first, second}));
}

/**
 * @brief Checks one resolution's factory table against the registry.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Every registration yields one entry with a name, a creator, a type
 *          key and a non-zero size; the table is the same object on every call.
 */
template <int W, int H> inline void verify_factory_table() {
  const std::vector<FactoryEntry> &table = hs_wasm::get_factory<W, H>();
  const std::vector<FactoryEntry> &again = hs_wasm::get_factory<W, H>();
  HS_EXPECT_EQ(table.size(), hs_wasm::EFFECT_REGISTRATIONS.size());
  HS_EXPECT_EQ(table.size(), static_cast<size_t>(HS_EFFECT_COUNT));
  HS_EXPECT_TRUE(&table == &again);

  for (size_t index = 0; index < table.size(); ++index) {
    const FactoryEntry &entry = table[index];
    HS_EXPECT_TRUE(!entry.identity.display_name.empty());
    HS_EXPECT_TRUE(!entry.identity.machine_id.empty());
    HS_EXPECT_TRUE(entry.identity.machine_id ==
                   hs_wasm::EFFECT_REGISTRATIONS[index].identity.machine_id);
    HS_EXPECT_TRUE(static_cast<bool>(entry.creator));
    HS_EXPECT_TRUE(entry.type_key != nullptr);
    HS_EXPECT_TRUE(entry.size > 0);
  }
}

/**
 * @brief Checks the name lookup setEffect() validates a UI string with.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Every HS_EFFECT_LIST name resolves to the entry carrying it; an
 *          unregistered name returns null.
 */
template <int W, int H> inline void verify_factory_lookup() {
  const auto lookup = [](std::string_view name) {
    return hs_wasm::find_factory_entry<W, H>(name);
  };
  for (const auto &registration : hs_wasm::EFFECT_REGISTRATIONS) {
    const char *name = registration.identity.display_name.data();
    const FactoryEntry *entry = lookup(name);
    HS_EXPECT_TRUE(entry != nullptr);
    if (entry)
      HS_EXPECT_TRUE(entry->identity.display_name == name);
  }
  for (const FactoryEntry &expected : hs_wasm::get_factory<W, H>()) {
    for (std::string_view key : expected.identity.accepted_ids())
      HS_EXPECT_TRUE(lookup(key) == &expected);
  }
  HS_EXPECT_TRUE(lookup("") == nullptr);
  HS_EXPECT_TRUE(lookup("NoSuchEffect") == nullptr);
  // A prefix of a registered name must miss: the scan compares whole names.
  const std::string_view first =
      hs_wasm::EFFECT_REGISTRATIONS[0].identity.display_name;
  HS_EXPECT_TRUE(lookup(first.substr(0, first.size() - 1)) == nullptr);
}

/**
 * @brief Builds, renders and destroys every registered effect at one resolution.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each instance is destroyed before the next is built: Effect permits
 *          one live instance.
 */
template <int W, int H> inline void drive_factory_lifecycle() {
  std::printf("  -- factory create/render/destroy @ %dx%d --\n", W, H);
  for (const FactoryEntry &entry : hs_wasm::get_factory<W, H>()) {
    effects_tests::reset_effect_globals();

    pin_frame_clock(0);
    std::unique_ptr<Effect> effect = entry.creator();
    HS_EXPECT_TRUE(effect != nullptr);
    if (!effect)
      continue;
    HS_EXPECT_EQ(effect->width(), W);
    HS_EXPECT_EQ(effect->height(), H);

    effect->init();
    HS_EXPECT_EQ(entry.preset_count, effect->getPresetCount());
    bool lit = false;
    for (int f = 0; f < FACTORY_OUTPUT_WINDOW; ++f) {
      pin_frame_clock(f);
      effect->draw_frame();
      // Consume the queued frame, else the next Canvas ctor spin-waits.
      effect->advance_display();
      for (int y = 0; y < H && !lit; ++y)
        for (int x = 0; x < W; ++x) {
          const Pixel p = effect->get_pixel(x, y);
          lit = lit || p.r || p.g || p.b;
        }
      if (lit && f + 1 >= FACTORY_FRAMES)
        break;
    }
    if (!lit)
      std::printf("  factory produced no output: %s @ %dx%d\n",
                  entry.identity.display_name.data(), W, H);
    HS_EXPECT_TRUE(lit);
    effect.reset();
    hs::clear_mock_time();
  }
}

/**
 * @brief Exercises the HS_RESOLUTIONS runtime dispatch.
 * @details Every declared row must dispatch and carry its own (W,H) into the
 *          callable; a size off the list must dispatch nowhere rather than fall
 *          into the nearest row.
 */
inline void test_resolution_dispatch() {
  for (const hs_wasm::WasmResolution &row : hs_wasm::WASM_RESOLUTIONS) {
    HS_EXPECT_TRUE(hs_wasm::wasm_resolution_supported(row.w, row.h));
    int seen_w = 0, seen_h = 0, hits = 0;
    const bool dispatched =
        hs_wasm::dispatch_resolution(row.w, row.h, [&]<int W, int H>() {
          seen_w = W;
          seen_h = H;
          ++hits;
        });
    HS_EXPECT_TRUE(dispatched);
    HS_EXPECT_EQ(hits, 1);
    HS_EXPECT_EQ(seen_w, row.w);
    HS_EXPECT_EQ(seen_h, row.h);
  }

  HS_EXPECT_FALSE(hs_wasm::wasm_resolution_supported(288.7, 144.2));
  HS_EXPECT_FALSE(hs_wasm::wasm_resolution_supported(288, 144.2));
  HS_EXPECT_FALSE(hs_wasm::wasm_resolution_supported(
      std::numeric_limits<double>::quiet_NaN(), 144));
  HS_EXPECT_FALSE(hs_wasm::wasm_resolution_supported(4294967584.0, 144));
  HS_EXPECT_TRUE(!hs_wasm::wasm_resolution_supported(0, 0));
  HS_EXPECT_TRUE(!hs_wasm::wasm_resolution_supported(-1, -1));
  HS_EXPECT_TRUE(!hs_wasm::wasm_resolution_supported(
      hs_wasm::WASM_RESOLUTIONS[0].w, hs_wasm::WASM_RESOLUTIONS[0].h + 1));
  int spurious = 0;
  const bool off_list =
      hs_wasm::dispatch_resolution(97, 21, [&]<int W, int H>() { ++spurious; });
  HS_EXPECT_TRUE(!off_list);
  HS_EXPECT_EQ(spurious, 0);
}

/**
 * @brief Checks the constructor's bootstrap name resolves at its bootstrap row.
 */
inline void test_bootstrap_rows() {
  HS_EXPECT_TRUE(hs_wasm::dispatch_resolution(
      hs_wasm::WASM_RESOLUTIONS[0].w, hs_wasm::WASM_RESOLUTIONS[0].h,
      []<int W, int H>() {
        HS_EXPECT_TRUE(
            (hs_wasm::find_factory_entry<W, H>(
                 hs_wasm::EFFECT_REGISTRATIONS[0].identity.display_name) !=
             nullptr));
      }));
}

/** @brief Checks registry-owned stable preset identities without Effect vtable cost. */
inline void test_fixed_preset_ids() {
  const FactoryEntry *lattice =
      hs_wasm::find_factory_entry<96, 20>("lattice-melt");
  HS_EXPECT_TRUE(lattice != nullptr);
  HS_EXPECT_TRUE(lattice && lattice->preset_id != nullptr);
  if (lattice && lattice->preset_id) {
    HS_EXPECT_EQ(lattice->preset_count, 2u);
    HS_EXPECT_TRUE(lattice->preset_id(0) == "open-curl");
    HS_EXPECT_TRUE(lattice->preset_id(1) == "dense-curl");
    HS_EXPECT_TRUE(lattice->preset_id(2).empty());
  }

  const FactoryEntry *mobius =
      hs_wasm::find_factory_entry<96, 20>("mobius-grid");
  HS_EXPECT_TRUE(mobius != nullptr);
  HS_EXPECT_TRUE(mobius && mobius->preset_id != nullptr);
  if (mobius && mobius->preset_id) {
    HS_EXPECT_EQ(mobius->preset_count, 2u);
    HS_EXPECT_TRUE(mobius->preset_id(0) == "mobius-grid");
    HS_EXPECT_TRUE(mobius->preset_id(1) == "mobius-grid-2");
  }

  const FactoryEntry *harmonics =
      hs_wasm::find_factory_entry<96, 20>("SphericalHarmonics");
  HS_EXPECT_TRUE(harmonics != nullptr);
  HS_EXPECT_TRUE(harmonics && harmonics->preset_id != nullptr);
  if (harmonics && harmonics->preset_id) {
    HS_EXPECT_EQ(harmonics->preset_count,
                 hs_test::effects_tests::SH_PRESET_MODES.size());
    if (harmonics->preset_count !=
        hs_test::effects_tests::SH_PRESET_MODES.size())
      return;
    for (size_t index = 0; index < harmonics->preset_count; ++index) {
      const auto [l, m] =
          SHMath::decode_lm(hs_test::effects_tests::SH_PRESET_MODES[index]);
      char expected[24];
      std::snprintf(expected, sizeof(expected), "sh-l%d-m%d", l, m);
      HS_EXPECT_TRUE(harmonics->preset_id(index) == expected);
    }
    HS_EXPECT_TRUE(
        harmonics->preset_id(hs_test::effects_tests::SH_PRESET_MODES.size())
            .empty());
  }

  const FactoryEntry *fishbowl =
      hs_wasm::find_factory_entry<96, 20>("Fishbowl");
  HS_EXPECT_TRUE(fishbowl != nullptr);
  HS_EXPECT_TRUE(fishbowl && fishbowl->preset_id != nullptr);
  if (fishbowl && fishbowl->preset_id) {
    HS_EXPECT_EQ(fishbowl->preset_count, 1u);
    HS_EXPECT_TRUE(fishbowl->preset_id(0) == "fire-trail");
    HS_EXPECT_TRUE(fishbowl->preset_id(1).empty());
  }
}

/**
 * @brief Checks every HS_RESOLUTIONS row's factory table and name lookup.
 */
inline void test_factory_tables() {
  for (const hs_wasm::WasmResolution &row : hs_wasm::WASM_RESOLUTIONS)
    HS_EXPECT_TRUE(
        hs_wasm::dispatch_resolution(row.w, row.h, []<int W, int H>() {
          verify_factory_table<W, H>();
          verify_factory_lookup<W, H>();
        }));
}

/** @brief Verifies Phantasm and the registry derive the same per-effect seed. */
inline void test_phantasm_seed_identity() {
  size_t index = 0;
#define HS_VERIFY_PHANTASM_SEED(name, duration_seconds)                        \
  do {                                                                         \
    const FactoryEntry *entry = hs_wasm::find_factory_entry<288, 144>(#name);  \
    HS_EXPECT_TRUE(entry != nullptr);                                          \
    const auto *golden = identity_golden(#name);                               \
    HS_EXPECT_TRUE(golden != nullptr);                                         \
    if (entry && golden)                                                       \
      HS_EXPECT_EQ((HS_PHANTASM_EFFECT_SEEDS<288, 144>[index]), golden->seed); \
    ++index;                                                                   \
  } while (false);
  HS_PHANTASM_EFFECT_LIST(HS_VERIFY_PHANTASM_SEED)
#undef HS_VERIFY_PHANTASM_SEED
}

/** @brief Each installation seed path restarts every roster effect's stream. */
inline void test_effect_stream_restarts() {
#define HS_VERIFY_EFFECT_STREAM(name)                                          \
  do {                                                                         \
    using E = name<96, 20>;                                                    \
    const auto *golden = identity_golden(#name);                               \
    HS_EXPECT_TRUE(golden != nullptr);                                         \
    if (!golden)                                                               \
      break;                                                                   \
    hs::Pcg32 expected(golden->seed);                                          \
    const uint32_t first = expected();                                         \
    hs::random()();                                                            \
    hs::seed_effect_stream<E>();                                               \
    HS_EXPECT_EQ(hs::random()(), first);                                       \
    hs::random()();                                                            \
    hs::seed_effect_stream(E::EFFECT_IDENTITY);                                \
    HS_EXPECT_EQ(hs::random()(), first);                                       \
    hs::random()();                                                            \
    hs::seed_effect_stream(golden->seed);                                      \
    HS_EXPECT_EQ(hs::random()(), first);                                       \
  } while (false);
  HS_EFFECT_LIST(HS_VERIFY_EFFECT_STREAM)
#undef HS_VERIFY_EFFECT_STREAM
}

/**
 * @brief Drives the create/render/destroy lifecycle at every HS_RESOLUTIONS
 *        row.
 * @details Rows above the small-aspect pixel count run only in the FULL tier.
 */
inline void test_factory_lifecycle() {
  constexpr int SMALL_PIXELS = effects_tests::SMALL_W * effects_tests::SMALL_H;
  const bool full = effects_tests::effects_full_suite();
  for (const hs_wasm::WasmResolution &row : hs_wasm::WASM_RESOLUTIONS) {
    if (row.w * row.h > SMALL_PIXELS && !full) {
      std::printf("  [TIER] factory %dx%d omitted; set HS_EFFECTS_FULL=1\n",
                  row.w, row.h);
      continue;
    }
    HS_EXPECT_TRUE(hs_wasm::dispatch_resolution(
        row.w, row.h, []<int W, int H>() { drive_factory_lifecycle<W, H>(); }));
  }
}

/**
 * @brief Module entry point: runs the factory and dispatch checks.
 * @return The module's failure count.
 */
inline int run_effect_factory_tests() {
  hs_test::ModuleFixture fixture("effect_factory");
  test_identity_contract();
  test_bootstrap_rows();
  test_resolution_dispatch();
  test_factory_tables();
  test_phantasm_seed_identity();
  test_effect_stream_restarts();
  test_fixed_preset_ids();
  test_factory_lifecycle();
  return fixture.result();
}

} // namespace effect_factory_tests
} // namespace hs_test
