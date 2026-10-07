"""PlatformIO pre-build hook: point sketch discovery at THIS env's .ino.

pioino.FindInoNodes globs only `$PROJECT_SRC_DIR/*.ino` and ignores
build_src_filter, so it finds nothing under targets/<X>/. This override returns
exactly the sketch for $PIOENV.
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
