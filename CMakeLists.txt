cmake_minimum_required(VERSION 3.21)
project(Holosphere CXX)

# Two build modes share this file, selected by the toolchain:
#   * Emscripten toolchain  -> WASM simulator module (daydream)   [EMSCRIPTEN]
#   * Native clang toolchain -> unit-test executable + CTest      [else]
# Use the presets in CMakePresets.json (cmake --preset wasm-release / tests).

# Enable ccache if available (both modes benefit).
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_FOUND})
  message(STATUS "ccache enabled: ${CCACHE_FOUND}")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(EMSCRIPTEN)
  # ===========================================================================
  # WASM simulator build (daydream)
  # ===========================================================================
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
  endif()

  set(TARGET_NAME holosphere_wasm)
  add_executable(${TARGET_NAME}
    targets/wasm/wasm.cpp
    core/memory.cpp
    core/reaction_graph.cpp)

  target_include_directories(${TARGET_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core)

  target_link_options(${TARGET_NAME} PRIVATE
    "-lembind"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sEXPORT_NAME=createHolosphereModule"
    "-sSTACK_SIZE=8192")

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WASM Debug build (-O0 -g -sASSERTIONS)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O0" "-g")
    target_link_options(${TARGET_NAME} PRIVATE "-O0" "-g" "-sASSERTIONS=1")
  else()
    message(STATUS "WASM Release build (-O3 -ffast-math -msimd128 -flto)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O3" "-ffast-math" "-msimd128" "-flto")
    target_link_options(${TARGET_NAME} PRIVATE "-O3" "-msimd128" "-flto")
  endif()

  # Emscripten emits a .js loader as the "executable".
  set_target_properties(${TARGET_NAME} PROPERTIES SUFFIX ".js")

  # Install the simulator module + a self-contained README into the daydream
  # checkout. Defaults to the sibling-checkout layout (../daydream); override
  # with -DDAYDREAM_DIR=/path/to/daydream. Destinations are absolute, so this
  # install is independent of CMAKE_INSTALL_PREFIX.
  #
  # If the directory does not exist (e.g. a standalone Holosphere clone), skip
  # the install rules with a notice instead of silently creating a stray
  # ../daydream. DAYDREAM_DIR is resolved at configure time — re-run cmake after
  # creating the checkout. (The `install` target still exists; it just no-ops.)
  set(DAYDREAM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../daydream"
      CACHE PATH "Destination daydream checkout for the installed WASM module")

  if(EXISTS "${DAYDREAM_DIR}")
    message(STATUS "Simulator module installs to DAYDREAM_DIR: ${DAYDREAM_DIR}")
    install(TARGETS ${TARGET_NAME} DESTINATION "${DAYDREAM_DIR}")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.wasm"
            DESTINATION "${DAYDREAM_DIR}")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
            DESTINATION "${DAYDREAM_DIR}")
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/docs/screenshots/"
            DESTINATION "${DAYDREAM_DIR}/docs/screenshots")
  else()
    message(STATUS
      "DAYDREAM_DIR '${DAYDREAM_DIR}' not found — skipping simulator-module "
      "install (set -DDAYDREAM_DIR=<path> and re-run cmake to enable).")
  endif()

else()
  # ===========================================================================
  # Native unit tests (Clang) + CTest
  # ===========================================================================
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
  endif()
  enable_testing()
  add_subdirectory(tests)

  # Opt-out local test gate: point git at .githooks so the pre-commit hook
  # (build + run the suite before each commit) is ON BY DEFAULT once the tests
  # are configured. Per-commit opt-out: HS_SKIP_TESTS=1 / --no-verify.
  # Disable the auto-enable entirely with -DHS_INSTALL_GIT_HOOKS=OFF.
  option(HS_INSTALL_GIT_HOOKS "Auto-enable the pre-commit test hook (core.hooksPath)" ON)
  if(HS_INSTALL_GIT_HOOKS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git"
     AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.githooks/pre-commit")
    execute_process(
      COMMAND git -C "${CMAKE_CURRENT_SOURCE_DIR}" config --local core.hooksPath .githooks
      RESULT_VARIABLE _hs_hook_rc OUTPUT_QUIET ERROR_QUIET)
    if(_hs_hook_rc EQUAL 0)
      message(STATUS "Pre-commit test hook enabled (core.hooksPath=.githooks). "
                     "Opt out: HS_SKIP_TESTS=1 / --no-verify; disable: -DHS_INSTALL_GIT_HOOKS=OFF.")
    endif()
  endif()
endif()
