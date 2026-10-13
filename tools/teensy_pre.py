"""PlatformIO pre-build hook: point sketch discovery at THIS env's .ino.

pioino.FindInoNodes globs only `$PROJECT_SRC_DIR/*.ino` and ignores
build_src_filter, so it finds nothing under targets/<X>/. This override returns
exactly the sketch for $PIOENV. The hook moves -fno-use-cxa-atexit from
BUILD_FLAGS to CXXFLAGS so C translation units do not receive it.
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
