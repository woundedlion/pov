"""PlatformIO pre-build hook: point sketch discovery at THIS env's .ino.

PlatformIO discovers the Arduino sketch ONLY by globbing
`$PROJECT_SRC_DIR/*.ino` at the top level (pioino.FindInoNodes) and IGNORES
build_src_filter, so with src_dir = repo root and the sketches under
targets/<X>/ it finds nothing.

This overrides FindInoNodes to return exactly this env's sketch (keyed on
$PIOENV); PlatformIO converts it to targets/<X>/<X>.ino.cpp for
build_src_filter. Only ONE sketch's setup()/loop() is ever compiled.
"""

import os

Import("env")

SKETCH = {
    "holosphere": os.path.join("targets", "Holosphere", "Holosphere.ino"),
    "holosphere_dma": os.path.join("targets", "Holosphere", "Holosphere.ino"),
    "phantasm": os.path.join("targets", "Phantasm", "Phantasm.ino"),
    "phantasm8": os.path.join("targets", "Phantasm", "Phantasm.ino"),
    "profile": os.path.join("targets", "Profile", "Profile.ino"),
    "profile_o3": os.path.join("targets", "Profile", "Profile.ino"),
    "bench": os.path.join("targets", "Bench", "Bench.ino"),
}

pioenv = env["PIOENV"]
if pioenv not in SKETCH:
    raise SystemExit(f"teensy_pre: no sketch mapping for env '{pioenv}'")
sketch_path = os.path.join(env["PROJECT_DIR"], SKETCH[pioenv])


def _find_ino(env):
    return [env.File(sketch_path)]


env.AddMethod(_find_ino, "FindInoNodes")
