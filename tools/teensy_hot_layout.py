"""PlatformIO pre-build hook: two-pass link that lays out hot flash code.

For an env linked by tools/phantasm.ld, the link becomes: a probe link with an
empty tools/hot_layout.py fragment, the layout plan written from the probe, and
the final link with that fragment. The probe never reads a previous plan, so
one source tree always yields one layout.
"""

import os
import sys

Import("env")

_ROOT = env.subst("$PROJECT_DIR")
sys.path.insert(0, os.path.join(_ROOT, "tools"))
import hot_layout  # noqa: E402

if os.path.basename(env.subst("$LDSCRIPT_PATH")) == "phantasm.ld":
    _DIR = os.path.join(env.subst("$BUILD_DIR"), "hot_layout")
    _PROBE_DIR = os.path.join(_DIR, "probe")
    _FINAL_DIR = os.path.join(_DIR, "final")
    env.Replace(
        HS_HOT_LAYOUT_PROBE_DIR=_PROBE_DIR,
        HS_HOT_LAYOUT_FINAL_DIR=_FINAL_DIR,
        HS_HOT_LAYOUT_PROBE_ELF=os.path.join(_PROBE_DIR, "probe.elf"),
        HS_HOT_LAYOUT_PROBE_MAP=os.path.join(_PROBE_DIR, "probe.map"))

    def _toolchain(env):
        return hot_layout.Toolchain.beside(env.subst("$CC"), env["ENV"].get("PATH"))

    def _write_probe_fragment(target, source, env):
        for directory in (_PROBE_DIR, _FINAL_DIR):
            os.makedirs(directory, exist_ok=True)
        with open(os.path.join(_PROBE_DIR, hot_layout.FRAGMENT_NAME), "w",
                  encoding="utf-8") as fragment:
            fragment.write("/* Probe link: no hot code layout. */\n")
        return 0

    def _plan(target, source, env):
        try:
            hot_layout.write_plan(
                _toolchain(env), env.subst("$HS_HOT_LAYOUT_PROBE_ELF"),
                env.subst("$HS_HOT_LAYOUT_PROBE_MAP"), hot_layout.roster(_ROOT),
                os.path.join(_FINAL_DIR, hot_layout.FRAGMENT_NAME), _ROOT)
        except hot_layout.LayoutError as exc:
            print(f"::error::hot-layout: {exc}")
            return 2
        return 0

    _LINK_ARGS = "$LINKFLAGS {extra} $__RPATH $SOURCES $_LIBDIRFLAGS $_LIBFLAGS"
    env.Replace(LINKCOM=[
        env.Action(_write_probe_fragment, None),
        "${TEMPFILE('$LINK -o $HS_HOT_LAYOUT_PROBE_ELF "
        + _LINK_ARGS.format(extra="-Wl,-Map,$HS_HOT_LAYOUT_PROBE_MAP "
                                  "-L$HS_HOT_LAYOUT_PROBE_DIR")
        + "','Linking hot-layout probe')}",
        env.Action(_plan, "Planning hot code layout"),
        "${TEMPFILE('$LINK -o $TARGET "
        + _LINK_ARGS.format(extra="-L$HS_HOT_LAYOUT_FINAL_DIR")
        + "','$LINKCOMSTR')}",
    ])
    env.Depends("$BUILD_DIR/${PROGNAME}.elf",
                [os.path.join(_ROOT, "tools", name)
                 for name in ("hot_layout.py", "teensy_hot_layout.py")])
