Import("env")
# ============================================================================
# newlib-nano for the Phantasm size build (docs/teensy_ci_gate_spec.md §4.1
# relief valve). nano swaps in the reduced libc/libstdc++ (libc_nano/
# libstdc++_nano): smaller stdio, malloc, and an integer-only printf whose float
# path is gated behind a `_printf_float` weak reference we deliberately do NOT
# request — the device formats no floats — so newlib's _dtoa_r + the %f/%g bignum
# helpers never link.
#
# --specs=nano.specs must reach the LINK step (that is where the library is
# selected); a flag placed only in build_flags reaches the compiler but not the
# linker, so nano silently does not engage. Add it to CCFLAGS (shared by C and
# C++ — a C++ command is `$CXXFLAGS $CCFLAGS`, so adding to CXXFLAGS too would
# include nano.specs twice on one command and gcc fatals: "spec 'link' already
# defined as nano_link") and to LINKFLAGS, each exactly once.
#
# Safe to mix here because PlatformIO compiles the Teensy core and FastLED from
# source under these same flags, so the whole image shares nano's _reent/stdio
# ABI — the usual "prebuilt core built against full newlib" hazard does not apply.
# ============================================================================
env.Append(CCFLAGS=["--specs=nano.specs"], LINKFLAGS=["--specs=nano.specs"])
