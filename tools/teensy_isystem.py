"""PlatformIO post: extra script: configure dependency warnings before compilation.

First-party code (core/ effects/ hardware/ targets/) keeps its -Wall -Wextra
warnings visible; the vendored FastLED (.pio/libdeps/<env>/) and Teensy core +
bundled libraries (PlatformIO packages dir) do not.

Project sources (projenv / env): third-party include dirs move from -I to
-isystem. The dir MUST be removed from -I, not merely also given as -isystem: a
relative -I and an absolute -isystem are distinct paths to gcc, and the -I
(searched first) would win and keep warning.

LDF library builders (FastLED, SPI): -w suppresses diagnostics in their own
.c/.cpp bodies. FrameworkArduino is built from the main environment and does
not receive this -w; its body warnings remain in the log.
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


# Project sources: demote third-party headers to -isystem (keep first-party -I).
if not sum(_demote_includes(build_env) for build_env in (projenv, env)):
    raise SystemExit(
        "teensy_isystem: demoted 0 third-party include dirs — no CPPPATH entry "
        "matched " + ", ".join(_THIRD_PARTY_ROOTS) + "; no include path lies under "
        "PlatformIO's packages or libdeps roots.")

# Library builders: their own source is third-party; disable its warnings.
for lib_builder in env.GetLibBuilders():
    lib_builder.env.Append(CCFLAGS=["-w"])
