/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file hash_constants.h
 * @brief Constants of the 32-bit FNV-1a hash. */
#include <cstdint>

namespace hs {
/// Offset basis of 32-bit FNV-1a.
inline constexpr uint32_t FNV1A_32_BASIS = 2166136261u;
/// Prime multiplier of 32-bit FNV-1a.
inline constexpr uint32_t FNV1A_32_PRIME = 16777619u;
} // namespace hs
