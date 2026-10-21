"""PlatformIO pre: extra script: run C, C++ and .S compiles through ccache.

Applies only when ccache is on the build PATH and the target is not
`compiledb`; otherwise the commands are unchanged. Links, archives and the
sketch's .ino preprocessing still call the toolchain directly.
CCACHE_DISABLE=1 bypasses the cache for one build.

The launcher sits inside `$( $)`, outside SCons build signatures, so enabling
or removing ccache reuses existing objects and object-cache entries.
"""

import shutil

Import("env")

COMPILE_COMMANDS = ("CCCOM", "CXXCOM", "ASPPCOM")

if ("compiledb" not in COMMAND_LINE_TARGETS
        and shutil.which("ccache", path=env["ENV"].get("PATH"))):
    env.Replace(**{com: "$( ccache $) " + env[com] for com in COMPILE_COMMANDS})
