"""PlatformIO pre-build hook: point sketch discovery at THIS env's .ino.

pioino.FindInoNodes globs only `$PROJECT_SRC_DIR/*.ino` and ignores
build_src_filter, so it finds nothing under targets/<X>/. This override returns
exactly the sketch for $PIOENV. The converted sketch's line markers name the
sketch by its project-relative path and drop the converter's temporary file, so
the .ino.cpp is identical across builds and checkouts. The hook moves
-fno-use-cxa-atexit from BUILD_FLAGS to CXXFLAGS so C translation units do not
receive it.
"""

import os
import re

Import("env")

build_flags = env.get("BUILD_FLAGS", [])
if isinstance(build_flags, str):
    build_flags = [build_flags]
env.Replace(BUILD_FLAGS=[re.sub(r"(?<!\S)-fno-use-cxa-atexit(?!\S)", "", flags)
                         for flags in build_flags])
env.Append(CXXFLAGS=["-fno-use-cxa-atexit"])

SKETCH = {
    "holosphere": os.path.join("targets", "Holosphere", "Holosphere.ino"),
    "holosphere_dma": os.path.join("targets", "Holosphere", "Holosphere.ino"),
    "phantasm": os.path.join("targets", "Phantasm", "Phantasm.ino"),
    "phantasm_rev12": os.path.join("targets", "Phantasm", "Phantasm.ino"),
    "phantasm8": os.path.join("targets", "Phantasm", "Phantasm.ino"),
    "profile": os.path.join("targets", "Profile", "Profile.ino"),
    "profile_o3": os.path.join("targets", "Profile", "Profile.ino"),
    "bench": os.path.join("targets", "Bench", "Bench.ino"),
    "bench_standalone": os.path.join("targets", "Bench", "Bench.ino"),
    "bench_master": os.path.join("targets", "Bench", "Bench.ino"),
}

pioenv = env["PIOENV"]
if pioenv not in SKETCH:
    raise SystemExit(f"teensy_pre: no sketch mapping for env '{pioenv}'")
sketch_path = os.path.join(env["PROJECT_DIR"], SKETCH[pioenv])


def _find_ino(env):
    return [env.File(sketch_path)]


env.AddMethod(_find_ino, "FindInoNodes")


# `# 12 "path" 1` from the preprocessor, `#line 12 "path"` from the converter.
LINE_MARKER = re.compile(r'^#(?:line)? \d+ "((?:[^"\\]|\\.)*)"')


def portable_line_markers(contents, sketch, project_dir):
    """Rename `sketch` to its project-relative path in line markers and drop
    markers that name any other file."""
    relative = os.path.relpath(sketch, project_dir).replace(os.sep, "/")
    target = os.path.normcase(os.path.normpath(sketch))
    lines = []
    for line in contents.splitlines(keepends=True):
        marker = LINE_MARKER.match(line)
        if marker:
            named = re.sub(r"\\(.)", r"\1", marker.group(1))
            if os.path.normcase(os.path.normpath(named)) != target:
                continue
            line = line[:marker.start(1)] + relative + line[marker.end(1):]
        lines.append(line)
    return "".join(lines)


convert_ino_to_cpp = env.ConvertInoToCpp


def _convert_ino(env):
    convert_ino_to_cpp()
    cpp_path = sketch_path + ".cpp"
    with open(cpp_path, encoding="utf-8", newline="") as cpp:
        contents = cpp.read()
    with open(cpp_path, "w", encoding="utf-8", newline="") as cpp:
        cpp.write(portable_line_markers(contents, sketch_path, env["PROJECT_DIR"]))


env.AddMethod(_convert_ino, "ConvertInoToCpp")
