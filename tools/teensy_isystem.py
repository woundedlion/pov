"""PlatformIO post-build hook: silence vendored-dependency warnings.

The firmware build runs -Wall -Wextra (the warning-hygiene ratchet, spec §7.2).
First-party code (core/ effects/ hardware/ targets/) must keep its warnings
visible. The vendored dependencies are a different matter: FastLED (under
.pio/libdeps/<env>/) and the Teensy core + its bundled libraries (under the
PlatformIO packages dir) are PINNED and gitignored, so their ~400 warnings are
noise we cannot fix in place (an edit is lost on the next `pio run`). This hook
removes that noise on both paths a warning can reach the log, each surgical so no
first-party diagnostic is lost:

1. Project sources (projenv / env): first-party TUs #include FastLED and Teensy
   core headers. Move those third-party include dirs from -I to -isystem so gcc
   treats them as system headers and drops their -Wregister / -Wnarrowing /
   -Wdeprecated-copy, while -Wall/-Wextra still cover first-party headers. The dir
   MUST be removed from -I, not merely also given as -isystem: a relative -I and
   an absolute -isystem are distinct paths to gcc, and the -I (searched first)
   would win and keep warning.

2. Library builders (FastLED / FrameworkArduino / SPI / ... compiling their OWN
   .c/.cpp): -isystem cannot reach a warning in the body of the file being
   compiled, only in its headers. These are third-party libraries we do not gate,
   so disable their warnings wholesale with -w (which also covers the C-vs-C++
   "valid for C++ but not C" cc1 notes the shared build_flags trigger on the core
   .c files). First-party code is built by projenv, never a library builder, so
   its -Wall/-Wextra coverage is untouched.
"""

import os

Import("env", "projenv")  # noqa: F821  (SCons globals injected by PlatformIO)

# A path is third-party if it lives under PlatformIO's libdeps or packages trees.
# The repo's own include dirs (., core, effects, hardware) match none of these.
_THIRD_PARTY_MARKERS = (
    os.sep + ".pio" + os.sep,        # .pio/libdeps/<env>/FastLED/src
    os.sep + ".platformio" + os.sep, # ~/.platformio/packages/framework-arduinoteensy/...
    os.sep + "packages" + os.sep,    # PLATFORMIO_CORE_DIR override that relocates packages/
)


def _is_third_party(path):
    norm = os.sep + os.path.normpath(path).strip(os.sep) + os.sep
    return any(marker in norm for marker in _THIRD_PARTY_MARKERS)


def _demote_includes(build_env):
    """Move third-party include dirs from CPPPATH (-I) to -isystem so gcc treats
    their headers as system headers and stops warning, keeping first-party -I."""
    kept = []
    isystem = []
    for entry in build_env.get("CPPPATH", []):
        resolved = os.path.normpath(build_env.subst(str(entry)))
        if _is_third_party(resolved):
            isystem += ["-isystem", resolved]
        else:
            kept.append(entry)
    if isystem:
        build_env.Replace(CPPPATH=kept)
        build_env.Append(CCFLAGS=isystem)


# 1. Project sources: demote third-party headers to -isystem (keep first-party -I).
for build_env in (projenv, env):
    _demote_includes(build_env)

# 2. Library builders: their own source is third-party; disable its warnings.
for lib_builder in env.GetLibBuilders():
    lib_builder.env.Append(CCFLAGS=["-w"])
