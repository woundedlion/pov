Import("env")
# Emit a linker map (firmware.map) into this env's build dir; it does not
# change the ELF/hex.
map_path = env.subst("$BUILD_DIR") + "/firmware.map"
env.Replace(HS_FIRMWARE_MAP_PATH=map_path)
env.Append(LINKFLAGS=["-Wl,-Map," + map_path])
env.NoCache("$BUILD_DIR/firmware.elf")
