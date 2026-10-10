/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file attributes.h
 * @brief Code and data placement macros: selective -O3 regions, the cold/flash
 *        family, and per-variable PROGMEM sections.
 */

#ifndef ARDUINO
// Off-device the Arduino storage qualifiers expand to nothing.
#ifndef DMAMEM
#define DMAMEM
#endif

#ifndef PROGMEM
#define PROGMEM
#endif

#ifndef FLASHMEM
#define FLASHMEM
#endif

#ifndef FASTRUN
#define FASTRUN
#endif
#endif

/**
 * @def HS_O3_BEGIN
 * @brief Opens a region whose function definitions compile at -O3 on the -Os
 *        device image; no-op for every other build.
 * @details The optimize pragma resets flags to defaults, so fast-math is
 *          restated. no-unswitch-loops stops GCC duplicating per-pixel loop
 *          bodies into ITCM.
 */
/** @def HS_O3_END
 *  @brief Closes an `HS_O3_BEGIN` region. */
/** @def HS_O3_FN
 *  @brief Single-function form of `HS_O3_BEGIN`/`HS_O3_END`. */
#if defined(ARDUINO) && defined(__GNUC__) && !defined(__clang__) &&            \
    defined(__OPTIMIZE_SIZE__)
#define HS_O3_BEGIN                                                            \
  _Pragma("GCC push_options") _Pragma(                                         \
      "GCC optimize(\"O3\", \"fast-math\", \"no-finite-math-only\", \"no-unswitch-loops\")")
#define HS_O3_END _Pragma("GCC pop_options")
#define HS_O3_FN                                                               \
  __attribute__((optimize("O3", "fast-math", "no-finite-math-only",            \
                          "no-unswitch-loops")))
#else
#define HS_O3_BEGIN
#define HS_O3_END
#define HS_O3_FN
#endif

/** @def HS_COLD
 *  @brief Places a setup-only free function in FLASH. noclone blocks the
 *         .constprop/.isra clones, which drop the section attribute and land
 *         in ITCM. */
/** @def HS_FLASH_MEMBER
 *  @brief Inline/template member variant of `HS_COLD`; `cold` emits a unique
 *         .text.unlikely.* section that tools/phantasm.ld routes to FLASH. */
/** @def HS_FLASH_INLINE
 *  @brief `HS_FLASH_MEMBER` without noinline, so free functions stay
 *         inlinable. */
/** @def HS_HOT_FLASH_MEMBER
 *  @brief `HS_FLASH_MEMBER` via .text.hot.*, without marking the code cold. */
/** @def HS_COLD_MEMBER
 *  @brief Alias of `HS_FLASH_MEMBER` (selective -O3 plus `cold`), for inline,
 *         template and COMDAT functions, member or free; unlike `HS_COLD` it
 *         carries the optimization override. */
/** @def HS_NOINLINE_NOCLONE
 *  @brief noinline/noclone; no section or optimization change. */
#if defined(__GNUC__) && !defined(__clang__)
#define HS_COLD FLASHMEM __attribute__((noinline, noclone))
#define HS_FLASH_MEMBER HS_O3_FN __attribute__((cold, noinline, noclone))
#define HS_FLASH_INLINE HS_O3_FN __attribute__((cold, noclone))
#define HS_HOT_FLASH_MEMBER HS_O3_FN __attribute__((hot, noinline, noclone))
#define HS_COLD_MEMBER HS_FLASH_MEMBER
#define HS_NOINLINE_NOCLONE __attribute__((noinline, noclone))
#else
#define HS_COLD FLASHMEM
#define HS_FLASH_MEMBER
#define HS_FLASH_INLINE
#define HS_HOT_FLASH_MEMBER
#define HS_COLD_MEMBER
#define HS_NOINLINE_NOCLONE __attribute__((noinline))
#endif

/** @def HS_HOT_INLINE
 *  @brief always_inline in optimized builds. At -O0 inlined bodies keep their
 *         own stack slots, so the helper stays an ordinary call there. */
#ifdef __OPTIMIZE__
#define HS_HOT_INLINE __attribute__((always_inline))
#else
#define HS_HOT_INLINE
#endif

/**
 * @def HS_PROGMEM_UNIQUE(name)
 * @brief Flash placement for a table defined in a header, in its own
 *        `.progmem.<name>` section.
 * @details Use this, never PROGMEM, for any COMDAT (inline/template) table:
 *          the shared ".progmem" section groups a TU's tables into one COMDAT,
 *          and ld can discard the whole group, resolving the other tables to
 *          address 0 without a diagnostic. tools/phantasm.ld globs
 *          `.progmem*`.
 * @param name Table identifier, stringized into the section name.
 */
#ifdef ARDUINO
#define HS_PROGMEM_UNIQUE(name) __attribute__((section(".progmem." #name)))
#else
#define HS_PROGMEM_UNIQUE(name)
#endif
