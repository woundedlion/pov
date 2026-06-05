@echo off
REM Thin wrapper over the canonical CMake preset. Equivalent to:
REM   cmake --preset tests && cmake --build --preset tests && ctest --preset tests
REM Requires cmake, ninja, and the emsdk clang (via the EMSDK env var or a
REM sibling ../emsdk checkout) on the system. No Visual Studio prompt needed.
cmake --preset tests || exit /b 1
cmake --build --preset tests || exit /b 1
ctest --preset tests
