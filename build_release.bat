@echo off
REM Thin wrapper over the canonical CMake presets. Equivalent to:
REM   cmake --preset wasm-release && cmake --build --preset wasm-release-install
REM The WASM build needs the Emscripten toolchain; ensure EMSDK is set (run
REM emsdk_env once) — we fall back to a sibling emsdk checkout if it is not.
if "%EMSDK%"=="" if exist "%~dp0..\emsdk\emsdk_env.bat" call "%~dp0..\emsdk\emsdk_env.bat"
cmake --preset wasm-release || exit /b 1
cmake --build --preset wasm-release-install
