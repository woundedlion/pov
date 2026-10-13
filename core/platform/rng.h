/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file rng.h
 * @brief Deterministic PRNG (Pcg32), hs::rand_* helpers and epoch/effect seeds.
 */

#include "core/platform/effect_identity.h"
#include <algorithm>
#include <cstdint>
#include <string_view>
#include <type_traits>

namespace hs {
/**
 * @brief Small deterministic PRNG (PCG XSH-RR 64/32) — the process-wide RNG.
 * @details Models a UniformRandomBitGenerator. Consume it only through hs::
 *          helpers: a \<random\> algorithm or distribution draws an
 *          implementation-defined number of times. Reference: pcg32 by
 *          Melissa O'Neill.
 */
class Pcg32 {
public:
  using result_type = uint32_t; ///< Output word type.
  /** @return Smallest value operator() returns. */
  static constexpr result_type min() { return 0u; }
  /** @return Largest value operator() returns. */
  static constexpr result_type max() { return 0xFFFFFFFFu; }

  /**
   * @brief Constructs the generator in the state seed() sets for `seed`.
   * @param seed Seed value.
   */
  explicit Pcg32(uint64_t seed = 1337u) { this->seed(seed); }

  /**
   * @brief Re-initializes the generator to the deterministic state for `s`.
   * @param s Seed value (mirrors std::mt19937::seed).
   */
  void seed(uint64_t s) {
    state = 0u;
    inc = (STREAM_SEQ << 1u) | 1u; // stream id must be odd
    (*this)();
    state += s;
    (*this)();
  }

  /** @return The next uniformly distributed 32-bit output; advances state. */
  result_type operator()() {
    uint64_t old = state;
    state = old * 6364136223846793005ULL + inc;
    uint32_t xorshifted = static_cast<uint32_t>(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = static_cast<uint32_t>(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((0u - rot) & 31u));
  }

private:
  uint64_t state = 0u;
  uint64_t inc = 0u;
  static constexpr uint64_t STREAM_SEQ = 0x14057b7ef767814fULL;
};

/**
 * @brief Derives the shared-RNG seed for effect epoch @p epoch.
 * @param epoch Roster index of the effect being built (beacon-synchronized on the device).
 * @return 1337 when epoch == 0 — the identity seed the determinism contract
 *         pins; otherwise a splitmix64-mixed value of (1337, epoch).
 * @details Integer-only, so device and host agree bit-for-bit.
 */
constexpr uint64_t epoch_seed(uint32_t epoch) {
  if (epoch == 0)
    return 1337u;
  uint64_t z = 1337u + uint64_t{epoch} * 0x9E3779B97F4A7C15ULL;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

/**
 * @brief Derives a roster-position-independent seed from an immutable RNG key.
 * @param effect_id Immutable effect RNG key.
 * @return 64-bit seed hashed from `effect_id`.
 */
constexpr uint64_t stable_effect_seed(std::string_view effect_id) {
  uint64_t hash = 1469598103934665603ULL;
  for (char value : effect_id) {
    hash ^= static_cast<uint8_t>(value);
    hash *= 1099511628211ULL;
  }
  uint64_t z = hash ^ 1337u;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

/**
 * @brief Returns the global deterministic random number generator.
 * @return Reference to the process-wide Pcg32, initially seeded with 1337.
 * @details Bit-identical between device and simulator, unlike the FastLED
 *          `random8()`/`random16()` and Arduino `random()` paths. Main-loop-only
 *          (non-atomic function-local static).
 */
inline Pcg32 &random() {
  static Pcg32 gen(1337);
  return gen;
}

/**
 * @brief Resets the effect RNG stream before construction.
 * @param seed Stable effect seed or synchronized epoch fallback.
 */
inline void seed_effect_stream(uint64_t seed) { random().seed(seed); }

/**
 * @brief Resets the effect RNG stream from its immutable compatibility key.
 * @param effect_id Immutable effect RNG key.
 */
inline void seed_effect_stream(std::string_view effect_id) {
  seed_effect_stream(stable_effect_seed(effect_id));
}

/**
 * @brief Resets an effect's RNG without changing its historical seed key.
 * @param identity Canonical effect and compatibility identity.
 */
inline void seed_effect_stream(const EffectIdentity &identity) {
  seed_effect_stream(identity.rng_seed_key);
}

/**
 * @brief Resets the RNG stream for a registered effect type.
 * @tparam EffectType Concrete effect with a declared effect identity.
 */
template <typename EffectType> inline void seed_effect_stream() {
  seed_effect_stream(effect_seed_key<EffectType>(EffectType::EFFECT_ID));
}

/**
 * @brief Maps a raw 32-bit RNG draw onto the half-open interval [0.0, 1.0).
 * @param value Raw 32-bit RNG draw.
 * @return A float in [0.0, 1.0), clamped just below 1.0f at the top band.
 * @details Top-band draws round to exactly 1.0f in float32; they clamp to the
 *          float just below 1.0f.
 */
inline float random_to_unit(uint32_t value) {
  float r = static_cast<float>(value) / static_cast<float>(UINT32_MAX);
  constexpr float JUST_BELOW_ONE = 0x1.fffffep-1f; // nextafterf(1.0f, 0.0f)
  return r > JUST_BELOW_ONE ? JUST_BELOW_ONE : r;
}

/**
 * @brief Generates a pseudo-random floating-point number in [0.0, 1.0).
 * @return A random float in the half-open range [0.0, 1.0).
 */
inline float rand_f() {
  using Rng = std::remove_reference_t<decltype(hs::random())>;
  static_assert(
      Rng::max() == 0xFFFFFFFFu,
      "rand_f()/random_to_unit assume a 32-bit-range RNG; update them "
      "if the global generator changes.");
  return random_to_unit(static_cast<uint32_t>(hs::random()()));
}

/**
 * @brief Generates a pseudo-random float in [min, max].
 * @param min Lower bound (inclusive).
 * @param max Upper bound (inclusive through rounding).
 * @return A random float in [min, max]; the top draws can round to max.
 */
inline float rand_f(float min, float max) {
  return min + rand_f() * (max - min);
}

/**
 * @brief Generates a pseudo-random integer within a specified range.
 * @param min The minimum value (inclusive).
 * @param max The maximum value (exclusive).
 * @return A random integer in [min, max), or min without consuming a random
 * draw when max <= min.
 * @note Uses `% span`, so the result is modulo-biased for ranges that do not
 * divide 2^32 evenly. The span is computed unsigned, so spans wider than
 * INT_MAX do not overflow.
 */
inline int rand_int(int min, int max) {
  if (max > min) {
    const uint32_t span =
        static_cast<uint32_t>(max) - static_cast<uint32_t>(min);
    return static_cast<int>(static_cast<uint32_t>(min) +
                            (hs::random()() % span));
  }
  return min;
}

/**
 * @brief Shuffles [first, last) using the process-wide RNG.
 * @tparam It Random-access iterator.
 * @param first Range begin.
 * @param last Range end.
 * @note Use instead of std::shuffle, whose permutation and draw count are
 * implementation-defined. Descending Fisher-Yates, exactly one
 * rand_int(0, i + 1) draw per step.
 */
template <typename It> inline void shuffle(It first, It last) {
  const int n = static_cast<int>(last - first);
  for (int i = n - 1; i > 0; --i) {
    const int j = rand_int(0, i + 1);
    std::iter_swap(first + i, first + j);
  }
}

} // namespace hs
