Import("env")
# Emit a linker map (firmware.map) into this env's build dir for size/layout
# analysis. The map does not change the ELF/hex. $BUILD_DIR resolves per-env.
env.Append(LINKFLAGS=["-Wl,-Map," + env.subst("$BUILD_DIR") + "/firmware.map"])
