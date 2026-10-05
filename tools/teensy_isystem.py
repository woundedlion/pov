"""PlatformIO post: extra script: configure dependency warnings before compilation.

The firmware build runs -Wall -Wextra.
First-party code (core/ effects/ hardware/ targets/) must keep its warnings
visible. The vendored dependencies are a different matter: FastLED (under
.pio/libdeps/<env>/) and the Teensy core + its bundled libraries (under the
PlatformIO packages dir) are pinned and gitignored. This hook demotes their
header diagnostics and suppresses LDF library-builder warnings while leaving
first-party diagnostics visible.

1. Project sources (projenv / env): first-party TUs #include FastLED and Teensy
   core headers. Move those third-party include dirs from -I to -isystem so gcc
   treats them as system headers and drops their -Wregister / -Wnarrowing /
   -Wdeprecated-copy, while -Wall/-Wextra still cover first-party headers. The dir
   MUST be removed from -I, not merely also given as -isystem: a relative -I and
   an absolute -isystem are distinct paths to gcc, and the -I (searched first)
   would win and keep warning.

2. LDF library builders (FastLED, SPI): -w suppresses diagnostics in their own
   .c/.cpp bodies. FrameworkArduino is built from the main environment and does
   not receive this -w. Its headers use -isystem, but its body warnings remain
   for the warning gate's third-party path filter. Core C command-line cc1 notes
   remain in the log; the gate's fileless-warning matcher selects cc1plus/ld.
"""

import os

Import("env", "projenv")

# A path is third-party if it lives under PlatformIO's libdeps or packages trees.
# The repo's own include dirs (., core, effects, hardware) match none of these.
_THIRD_PARTY_ROOTS = tuple(
    os.path.normcase(os.path.realpath(env.subst(variable)))
    for variable in ("$PROJECT_PACKAGES_DIR", "$PROJECT_LIBDEPS_DIR")
)


def _is_third_party(path):
    norm = os.path.normcase(os.path.realpath(path))
    for root in _THIRD_PARTY_ROOTS:
        try:
            if os.path.commonpath((norm, root)) == root:
                return True
        except ValueError:
            pass
    return False


def _demote_includes(build_env):
    """Move third-party include dirs from CPPPATH (-I) to -isystem so gcc treats
    their headers as system headers and stops warning, keeping first-party -I.

    Returns the number of dirs demoted."""
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
    return len(isystem) // 2


# 1. Project sources: demote third-party headers to -isystem (keep first-party -I).
# Every Teensy build carries the framework core dir (…/packages/framework-
# arduinoteensy/cores/teensy4) on CPPPATH, so demoting nothing means the marker
# set stopped matching PlatformIO's layout — not that the paths are clean. Fail
# instead of silently reverting to a warning-flooded log: the warning gate's
# first-party filter drops the flood, so nothing else would notice.
if not sum(_demote_includes(build_env) for build_env in (projenv, env)):
    raise SystemExit(
        "teensy_isystem: demoted 0 third-party include dirs — no CPPPATH entry "
        "matched " + ", ".join(_THIRD_PARTY_ROOTS) + "; the vendored-path "
        "markers no longer match PlatformIO's layout.")

# 2. Library builders: their own source is third-party; disable its warnings.
for lib_builder in env.GetLibBuilders():
    lib_builder.env.Append(CCFLAGS=["-w"])
