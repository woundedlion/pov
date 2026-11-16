/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file float_bits.h
 * @brief Bitwise float comparisons for composed-effect tests. */

#include <bit>
#include <cstdint>

namespace hs_test::composed_effect_tests {

/** @brief Returns the stored bit pattern of a float.
 * @param value Float value.
 * @return Unsigned bit representation. */
inline uint32_t bits(float value) { return std::bit_cast<uint32_t>(value); }

} // namespace hs_test::composed_effect_tests
