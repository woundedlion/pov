/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// No include guard: each inclusion expands the playlist once.

/**
 * @file phantasm_playlist_rows.h
 * @brief Expands every phantasm_playlist.def row through the includer's
 *        HS_PHANTASM_ROW and HS_PHANTASM_EXCLUDED_ROW, then undefines both.
 *
 * Define either macro before including; an undefined one expands to nothing.
 *   HS_PHANTASM_ROW(name, duration_seconds)  one per show row, in show order
 *   HS_PHANTASM_EXCLUDED_ROW(name)           one per exclude row
 * A show row without seconds passes hs_cycle_seconds<name>().
 */

#ifndef HS_PHANTASM_ROW
#define HS_PHANTASM_ROW(name, duration_seconds)
#endif
#ifndef HS_PHANTASM_EXCLUDED_ROW
#define HS_PHANTASM_EXCLUDED_ROW(name)
#endif

/**
 * @brief Selects the show-row expansion by argument count.
 * @param name Effect class name.
 * @param seconds Show duration, or the fixed-duration expansion when absent.
 * @param row The expansion that matches the argument count.
 */
#define HS_PHANTASM_SHOW_PICK(name, seconds, row, ...) row
/**
 * @brief Expands a show row with a fixed duration.
 * @param name Effect class name.
 * @param seconds Show duration in seconds.
 */
#define HS_PHANTASM_SHOW_FOR(name, seconds) HS_PHANTASM_ROW(name, seconds)
/**
 * @brief Expands a show row that plays one cycle of the effect.
 * @param name Effect class name.
 */
#define HS_PHANTASM_SHOW_CYCLE(name)                                           \
  HS_PHANTASM_ROW(name, hs_cycle_seconds<name>())

#pragma push_macro("show")
#pragma push_macro("exclude")
#undef show
#undef exclude
/** @brief A playlist show row: `show(Effect, seconds)` or `show(Effect)`. */
#define show(...)                                                              \
  HS_PHANTASM_SHOW_PICK(__VA_ARGS__, HS_PHANTASM_SHOW_FOR,                     \
                        HS_PHANTASM_SHOW_CYCLE, )(__VA_ARGS__)
/**
 * @brief A playlist exclude row.
 * @param name Effect class name.
 */
#define exclude(name) HS_PHANTASM_EXCLUDED_ROW(name)

#include "targets/Phantasm/phantasm_playlist.def"

#undef show
#undef exclude
#pragma pop_macro("exclude")
#pragma pop_macro("show")
#undef HS_PHANTASM_SHOW_CYCLE
#undef HS_PHANTASM_SHOW_FOR
#undef HS_PHANTASM_SHOW_PICK
#undef HS_PHANTASM_EXCLUDED_ROW
#undef HS_PHANTASM_ROW
