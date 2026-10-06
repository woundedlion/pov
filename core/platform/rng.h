/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file rng.h
 * @brief Deterministic PRNG (Pcg32), hs::rand_* helpers, epoch/effect seeds
 *        and the stable effect identity used by effect seeding.
 */

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
  using result_type = uint32_t;
  static constexpr result_type min() { return 0u; }
  static constexpr result_type max() { return 0xFFFFFFFFu; }

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
 * @brief Returns an effect type's persisted ID, falling back to its class name.
 * @tparam EffectType Concrete effect type at a fixed resolution.
 * @param class_name Effect class name.
 */
template <typename EffectType>
constexpr std::string_view stable_effect_id(std::string_view class_name) {
  if constexpr (requires { EffectType::EFFECT_ID; })
    return EffectType::EFFECT_ID;
  return class_name;
}

/** @brief Derives a roster-position-independent seed from an effect ID. */
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
 * @brief Maps a raw RNG draw in [0, max] onto the half-open interval [0.0, 1.0).
 * @param value Raw RNG draw, in [0, max].
 * @param max Maximum possible draw value.
 * @pre max == UINT32_MAX.
 * @return A float in [0.0, 1.0), clamped just below 1.0f at the top band.
 * @details Top-band draws round to exactly 1.0f in float32; they clamp to the
 *          float just below 1.0f.
 */
inline float random_to_unit(uint32_t value, uint32_t max) {
  float r = static_cast<float>(value) / static_cast<float>(max);
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
  return random_to_unit(static_cast<uint32_t>(hs::random()()),
                        static_cast<uint32_t>(Rng::max()));
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
