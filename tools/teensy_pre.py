"""PlatformIO pre-build hook: point sketch discovery at THIS env's .ino.

Phase-0 finding (corrects spec §5/§6): PlatformIO discovers the Arduino sketch
ONLY by globbing `$PROJECT_SRC_DIR/*.ino` at the top level (pioino.FindInoNodes)
and IGNORES build_src_filter. With src_dir = repo root and the two sketches under
targets/<X>/, that glob finds nothing, so setup()/loop() never link.

We keep src_dir = the repo root — so the shared core/engine/*.cpp build as NORMAL project
sources and inherit the LDF-resolved library include paths (FastLED, and the
framework's SPI) they need — and simply override FindInoNodes to return exactly
this env's sketch. PlatformIO then converts it to targets/<X>/<X>.ino.cpp, which
build_src_filter picks up (see platformio.ini). Selecting the sketch here (keyed
on $PIOENV) also guarantees only ONE sketch's setup()/loop() is ever compiled,
even though both .ino files define them.
"""

import os

Import("env")  # noqa: F821  (SCons global injected by PlatformIO)

SKETCH = {
    "holosphere": os.path.join("targets", "Holosphere", "Holosphere.ino"),
    "holosphere_dma": os.path.join("targets", "Holosphere", "Holosphere.ino"),
    "phantasm": os.path.join("targets", "Phantasm", "Phantasm.ino"),
    "phantasm8": os.path.join("targets", "Phantasm", "Phantasm.ino"),
    "profile": os.path.join("targets", "Profile", "Profile.ino"),
    "profile_o3": os.path.join("targets", "Profile", "Profile.ino"),
}

pioenv = env["PIOENV"]
if pioenv not in SKETCH:
    raise SystemExit(f"teensy_pre: no sketch mapping for env '{pioenv}'")
sketch_path = os.path.join(env["PROJECT_DIR"], SKETCH[pioenv])


def _find_ino(env):
    return [env.File(sketch_path)]


env.AddMethod(_find_ino, "FindInoNodes")
