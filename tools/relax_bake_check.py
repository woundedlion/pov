"""PlatformIO pre-build freshness gate for generated relax bakes."""

import os
import subprocess
import sys


Import("env")  # noqa: F821

toolchain_dir = env.PioPlatform().get_package_dir(
    "toolchain-gccarmnoneeabi-teensy"
)
compiler_name = "arm-none-eabi-g++.exe" if os.name == "nt" else "arm-none-eabi-g++"
compiler_path = os.path.join(toolchain_dir, "bin", compiler_name)
subprocess.run(
    [
        sys.executable,
        "tools/relax_bakes.py",
        "compiler",
        "--path",
        compiler_path,
    ],
    cwd=env["PROJECT_DIR"],
    check=True,
)

external_flags = os.environ.get("PLATFORMIO_BUILD_FLAGS", "")
build_flags = env.get("BUILD_FLAGS", [])
resolved_flags = (
    " ".join(str(flag) for flag in build_flags)
    if isinstance(build_flags, (list, tuple))
    else str(build_flags)
)
if "HS_RELAX_BAKE_EXTRACT" not in f"{external_flags} {resolved_flags}":
    subprocess.run(
        [sys.executable, "tools/relax_bakes.py", "check"],
        cwd=env["PROJECT_DIR"],
        check=True,
    )
