# Holosphere native unit tests.
#
# The quick canonical runner is ../build_tests.bat. This file additionally
# exposes the suite to CMake/CTest so CI can build and run it portably:
#
#   cmake -S tests -B build_test -G Ninja
#   cmake --build build_test
#   ctest --test-dir build_test --output-on-failure
#
# The root CMakeLists.txt (a WASM/Emscripten build) does add_subdirectory(tests);
# the unit tests are native-only, so we no-op under Emscripten to leave the WASM
# build untouched.
#
# NOTE: on Windows with a clang toolchain, CMake's platform module also wants a
# resource compiler — pass -DCMAKE_RC_COMPILER=<path-to-llvm-rc>. Linux/macOS CI
# need nothing extra. The verified quick runner on this repo is build_tests.bat.
if(EMSCRIPTEN)
  return()
endif()

cmake_minimum_required(VERSION 3.16)

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
