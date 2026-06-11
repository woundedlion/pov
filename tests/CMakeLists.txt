# Holosphere native unit tests.
#
# Canonical flow is the `tests` CMake preset (see CMakePresets.json), which uses
# cmake/toolchain-native-clang.cmake to locate Clang and, on Windows, the SDK
# rc.exe automatically — no Visual Studio Developer Prompt or manual flags:
#
#   cmake --preset tests
#   cmake --build --preset tests
#   ctest --preset tests
#
# build_tests.bat is a thin wrapper over exactly those three commands.
#
# The root CMakeLists.txt (a WASM/Emscripten build) does add_subdirectory(tests);
# the unit tests are native-only, so we no-op under Emscripten to leave the WASM
# build untouched. (A standalone `cmake -S tests -B build_test` configure also
# works for CI, but the preset is the supported path.)
if(EMSCRIPTEN)
  return()
endif()

# Matches the root floor: the native-clang toolchain sets CMAKE_LINKER_TYPE
# (added in CMake 3.29), which is silently ignored on older CMake. This applies
# to the standalone `cmake -S tests` path too, so gate it here as well.
cmake_minimum_required(VERSION 3.29)

# Allow standalone configuration (cmake -S tests) as well as being pulled in by
# a native root build.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(holosphere_tests CXX)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(run_tests
  ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/reaction_graph.cpp)

# HS_TEST_BUILD widens the test-only inline/arena budgets (MSVC STL's
# std::function is larger than libc++/inplace_function on the device); the
# device footprint is unchanged. See core/animation.h and core/memory.h.
target_compile_definitions(run_tests PRIVATE HS_TEST_BUILD)

# Headers are included as "tests/..." and "core/..." from the repo root.
target_include_directories(run_tests PRIVATE ${HS_ROOT} ${HS_ROOT}/core)

enable_testing()
add_test(NAME unit_tests COMMAND run_tests)
