"""PlatformIO post: extra script: silence third-party warnings before compilation.

Third-party include dirs (PlatformIO packages and libdeps) move from -I to
-isystem. Each MUST be removed from -I, not merely added as -isystem: gcc sees a
relative -I and an absolute -isystem as distinct paths, and the -I wins.

LDF library builders (FastLED, SPI) get -w. FrameworkArduino is built from the
main environment and keeps its body warnings.
"""

import os

Import("env", "projenv")

# A path is third-party if it lives under PlatformIO's libdeps or packages trees.
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
