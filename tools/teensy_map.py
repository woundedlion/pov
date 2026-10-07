Import("env")
# Emit a linker map (firmware.map) into this env's build dir; it does not
# change the ELF/hex.
env.Append(LINKFLAGS=["-Wl,-Map," + env.subst("$BUILD_DIR") + "/firmware.map"])
