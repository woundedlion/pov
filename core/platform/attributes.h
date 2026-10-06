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

// ---------------------------------------------------------------------------
// HS_O3_BEGIN / HS_O3_END: compile the enclosed function definitions at -O3 on
// the -Os device image; no-op for every other build. HS_O3_FN is the
// single-function form.
// The optimize pragma resets flags to defaults, so fast-math is restated.
// no-unswitch-loops stops GCC duplicating per-pixel loop bodies into ITCM.
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// HS_COLD: places a setup-only free function in FLASH. noclone blocks the
// .constprop/.isra clones, which drop the section attribute and land in ITCM.
// HS_FLASH_MEMBER: inline/template member variant; `cold` emits a unique
// .text.unlikely.* section that tools/phantasm.ld routes to FLASH.
// HS_HOT_FLASH_MEMBER: the same via .text.hot.*, without marking the code cold.
// HS_FLASH_INLINE omits noinline so free functions stay inlinable.
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// HS_PROGMEM_UNIQUE: flash placement for a table defined in a header. Use this,
// never PROGMEM, for any COMDAT (inline/template) table: the shared ".progmem"
// section groups a TU's tables into one COMDAT, and ld can discard the whole
// group, resolving the other tables to address 0 without a diagnostic.
// tools/phantasm.ld globs `.progmem*`.
// ---------------------------------------------------------------------------
#ifdef ARDUINO
#define HS_PROGMEM_UNIQUE(name) __attribute__((section(".progmem." #name)))
#else
#define HS_PROGMEM_UNIQUE(name)
#endif
