/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Buffer watchdog override parsing for native tests. */

#include <charconv>
#include <cstdint>
#include <cstring>

namespace hs_test {

/** @brief Reads the buffer watchdog override in microseconds.
 * @param text Unsigned decimal override, or null.
 * @return The override, or zero for an absent or invalid value.
 */
inline uint32_t parse_buffer_watchdog_us(const char *text) {
  if (!text)
    return 0;
  uint32_t value = 0;
  const char *end = text + std::strlen(text);
  const auto parsed = std::from_chars(text, end, value);
  return parsed.ec == std::errc{} && parsed.ptr == end ? value : 0;
}

} // namespace hs_test
