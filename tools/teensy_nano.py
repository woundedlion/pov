Import("env")
# ============================================================================
# newlib-nano for the Phantasm size build: the reduced libc/libstdc++ with an
# integer-only printf (no `_printf_float` reference is requested, so newlib's
# float formatting never links).
#
# --specs=nano.specs must reach the LINK step; a flag only in build_flags never
# reaches the linker. It goes to CCFLAGS (a C++ command is
# `$CXXFLAGS $CCFLAGS`, so adding it to CXXFLAGS too repeats it and gcc fatals:
# "spec 'link' already defined as nano_link") and to LINKFLAGS, each exactly
# once.
#
# The Teensy core and FastLED compile from source under these same flags, so the
# whole image shares nano's _reent/stdio ABI.
# ============================================================================
env.Append(CCFLAGS=["--specs=nano.specs"], LINKFLAGS=["--specs=nano.specs"])
