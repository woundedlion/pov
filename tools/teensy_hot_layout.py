"""PlatformIO pre-build hook: two-pass link that lays out hot flash code.

For an env linked by tools/phantasm.ld, the link becomes: a probe link with an
empty tools/hot_layout.py fragment, the layout plan written from the probe, the
final link with that fragment, and the hot-layout gate on the final image. The
probe never reads a previous plan, so one source tree always yields one layout.
"""

import os
import re
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

    def _final_map(env):
        for flag in env.subst("$LINKFLAGS").split():
            if flag.startswith("-Wl,-Map,"):
                return flag[len("-Wl,-Map,"):]
        raise hot_layout.LayoutError("the final link writes no -Wl,-Map")

    def _sketch(source):
        for node in source:
            match = re.search(r"([A-Za-z]+)\.ino\.cpp\.o$", str(node))
            if match:
                return match.group(1)
        return ""

    def _expected(env, source):
        sketch = _sketch(source)
        if sketch == "Phantasm":
            return hot_layout.phantasm_playlist(_ROOT)
        if sketch == "Profile":
            for define in env.get("CPPDEFINES", []):
                if isinstance(define, (list, tuple)) and define[0] == "HS_PROFILE_TARGET":
                    return [str(define[1])]
            with open(os.path.join(_ROOT, "targets/Profile/Profile.ino"),
                      encoding="utf-8") as sketch_file:
                text = sketch_file.read()
            return re.findall(r"#define HS_PROFILE_TARGET (\w+)", text)[:1]
        return []

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

    def _gate(target, source, env):
        try:
            lines, passed = hot_layout.check(
                _toolchain(env), str(target[0]), _final_map(env),
                hot_layout.roster(_ROOT), _expected(env, source), env["PIOENV"])
        except hot_layout.LayoutError as exc:
            print(f"::error::hot-layout: {exc}")
            return 2
        with open(os.path.join(_DIR, "report.txt"), "w", encoding="utf-8") as report:
            report.write("\n".join(lines) + "\n")
        print("\n".join(lines))
        if not passed:
            print(f"::error::hot-layout: {env['PIOENV']} failed the hot-layout gate; "
                  f"see {_DIR}/report.txt")
            return 1
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
        env.Action(_gate, "Checking hot code layout"),
    ])
    env.Depends("$BUILD_DIR/firmware.elf",
                [os.path.join(_ROOT, "tools", name)
                 for name in ("cpp_roster.py", "hot_layout.py", "teensy_gate.py",
                              "teensy_hot_layout.py", "phantasm.ld")])
