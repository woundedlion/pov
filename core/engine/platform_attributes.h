/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file platform_attributes.h
 * @brief Code and data placement macros: selective -O3 regions, the cold/flash
 *        family, and per-variable PROGMEM sections.
 */

#ifndef ARDUINO
// Off-device the Arduino storage qualifiers have no meaning; supply empty
// definitions so the placement macros below expand everywhere. On device they
// come from Arduino.h.
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
// the -Os device image (selective optimization).
// Active only for device GCC building at -Os (__OPTIMIZE_SIZE__): the holosphere
// -O3 image, host clang, and WASM see no-ops, so those builds are byte-identical.
// The fast-math flags are restated because GCC 11's optimize pragma rebuilds
// optimization flags from defaults, dropping the command-line ones for the
// region. no-unswitch-loops keeps GCC 15 from emitting a second copy of the
// region's per-pixel loop bodies, which only overflows ITCM: one copy ever runs.
// HS_O3_FN is the shared single-function attribute and backs HS_COLD_MEMBER.
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
// HS_COLD: keep a setup-only function off the fast ITCM banks. FLASHMEM routes it
// to FLASH; noinline collapses per-call-site inline copies and noclone blocks the
// .constprop/.isra IPA clones (which drop the section attribute and land in ITCM
// regardless). Apply ONLY to internal-linkage (`static`) free functions on cold
// paths (mesh/solid construction): a section attribute on a COMDAT (inline/template
// member) function is a section-type conflict. Off-device it degrades to a no-op.
// HS_FLASH_MEMBER is the COMDAT-safe variant for inline/template member
// functions: GCC's `cold` attribute supplies a unique .text.unlikely.* section,
// and tools/phantasm.ld routes that section to FLASH. HS_HOT_FLASH_MEMBER uses
// the corresponding .text.hot.* route for measured hot code that executes from
// cached flash without telling the optimizer it is cold. HS_COLD_MEMBER names the
// setup-only use; HS_FLASH_MEMBER also supports explicitly measured code
// placement. On the -Os device image both use HS_O3_FN. HS_FLASH_INLINE is the
// variant for a free function declared `inline`, which GCC's -Wattributes
// rejects the noinline on; `cold` alone still supplies the .text.unlikely.*
// section, and a variadic [[noreturn]] body is not an inline candidate anyway.
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
// never PROGMEM, for any COMDAT (inline/template) table.
//
// Teensy's PROGMEM is the single fixed section ".progmem", so every inline
// PROGMEM variable in a translation unit lands in one section inside one COMDAT
// group, signed by whichever symbol sits at offset 0. Two translation units that
// pick the same signature make ld discard a whole group: the other tables in it
// go too, their weak references resolve to address 0, and the first read faults
// with no linker diagnostic. A per-variable section name gives each table its
// own group, so tables dedupe individually. tools/phantasm.ld globs `.progmem*`.
// ---------------------------------------------------------------------------
#ifdef ARDUINO
#define HS_PROGMEM_UNIQUE(name) __attribute__((section(".progmem." #name)))
#else
#define HS_PROGMEM_UNIQUE(name)
#endif
