Import("env")
# newlib-nano (integer-only printf) for the Phantasm size build.
# --specs=nano.specs goes to CCFLAGS and LINKFLAGS exactly once each: a C++
# command is `$CXXFLAGS $CCFLAGS`, so adding it to CXXFLAGS too makes gcc fatal
# with "spec 'link' already defined as nano_link".
env.Append(CCFLAGS=["--specs=nano.specs"], LINKFLAGS=["--specs=nano.specs"])
