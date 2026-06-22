Import("env")
# Emit a linker map (firmware.map) into this env's build dir for size/layout
# analysis. The map is a separate text artifact under .pio/ (gitignored): it does
# NOT change the binary (ELF/hex are byte-identical with or without it), costs
# well under a second at link, and is the per-symbol/per-object source the size
# investigations read. $BUILD_DIR resolves per-env, so each target gets its own.
env.Append(LINKFLAGS=["-Wl,-Map," + env.subst("$BUILD_DIR") + "/firmware.map"])
