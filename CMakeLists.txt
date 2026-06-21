# 3.29 is the real floor: CMakePresets schema v6 needs >=3.25 (3.21-3.24 reject
# the file outright) and the native test toolchain's CMAKE_LINKER_TYPE was added
# in 3.29 (silently ignored on 3.25-3.28, so the Windows tests build falls back
# to the wrong linker — the exact failure toolchain-native-clang.cmake exists to
# prevent). Declare it honestly so an old CMake fails loudly instead of mis-linking.
cmake_minimum_required(VERSION 3.29)
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
    message(STATUS "WASM Debug build (-O0 -g -sASSERTIONS, 64 KB stack)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O0" "-g")
    # Override the shared 8 KB STACK_SIZE for debug only. -O0 disables inlining
    # and stack-slot coalescing, inflating frame sizes severalfold over the -O3
    # release build, so the tight release budget overflows under no optimization.
    # The debug preset is never exercised by CI or the smoke script, so size it
    # generously rather than rely on -sASSERTIONS to merely report the abort.
    # (Last -sSTACK_SIZE on the link line wins; release is unaffected.)
    target_link_options(${TARGET_NAME} PRIVATE "-O0" "-g" "-sASSERTIONS=1"
      "-sSTACK_SIZE=65536")
  else()
    # -fno-finite-math-only is required AFTER -ffast-math: -ffast-math implies
    # -ffinite-math-only, which folds std::isfinite() to constant true and lets
    # the compiler assume no NaN/Inf module-wide. The render sink relies on real
    # finite semantics (filter.h: a NaN coord must not silently UB-cast and
    # bypass the clip/wrap guards). We keep the rest of -ffast-math (reassoc,
    # FMA, reciprocal, -fno-math-errno) for the SIMD float throughput.
    message(STATUS "WASM Release build (-O3 -ffast-math -fno-finite-math-only -msimd128 -flto)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O3" "-ffast-math" "-fno-finite-math-only" "-msimd128" "-flto")
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

    # Record the engine commit this WASM module was built from, plus the module's
    # own content hash, alongside it. The daydream deploy gate reads the SHA to pin
    # its engine checkout to the exact source the committed module corresponds to
    # (instead of a decoupled engine HEAD), and asserts the committed .wasm matches
    # the recorded hash so a stale or swapped binary can't ride a green gate. Both
    # are captured at install time (immediately after the build) from the same
    # build event, and committed together with the .wasm they describe.
    install(CODE "
      execute_process(
        COMMAND git -C \"${CMAKE_CURRENT_SOURCE_DIR}\" rev-parse HEAD
        OUTPUT_VARIABLE HS_ENGINE_SHA
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE HS_GIT_RC)
      # Dirty-tree guard: a module built from uncommitted (tracked) source records
      # provenance the deploy gate cannot reproduce — the gate checks out the bare
      # SHA (clean) and tests that, while the shipped binary came from dirty source.
      # Mark the SHA -dirty (git describe convention) so the gate's checkout of that
      # ref fails hard instead of deploying an unverifiable artifact. Untracked files
      # don't count (build output, node_modules); only tracked modifications do.
      execute_process(
        COMMAND git -C \"${CMAKE_CURRENT_SOURCE_DIR}\" status --porcelain --untracked-files=no
        OUTPUT_VARIABLE HS_GIT_DIRTY
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      if(NOT \"\${HS_GIT_DIRTY}\" STREQUAL \"\")
        set(HS_ENGINE_SHA \"\${HS_ENGINE_SHA}-dirty\")
        message(WARNING \"Engine tree has uncommitted changes; recording ${TARGET_NAME}.sha as \${HS_ENGINE_SHA}. The deploy gate will reject this ref — rebuild from a clean (committed) tree before deploying.\")
      endif()
      if(HS_GIT_RC EQUAL 0)
        file(WRITE \"${DAYDREAM_DIR}/${TARGET_NAME}.sha\" \"\${HS_ENGINE_SHA}\\n\")
        message(STATUS \"Wrote engine provenance ${TARGET_NAME}.sha: \${HS_ENGINE_SHA}\")
      else()
        # Fail-fast at the violation site: a missing .sha but a written .wasm.sha256
        # is a half-written provenance trio that the deploy gate rejects anyway, so
        # abort the install here rather than returning success with an incomplete
        # set (the .wasm.sha256 write below is never reached).
        message(FATAL_ERROR \"git rev-parse failed (rc=\${HS_GIT_RC}); aborting install so ${TARGET_NAME} provenance is not left half-written (the deploy gate needs the .sha)\")
      endif()

      # sha256sum -c format ('<hash>  <filename>'), checked from the daydream root.
      file(SHA256 \"${DAYDREAM_DIR}/${TARGET_NAME}.wasm\" HS_WASM_SHA256)
      file(WRITE \"${DAYDREAM_DIR}/${TARGET_NAME}.wasm.sha256\"
           \"\${HS_WASM_SHA256}  ${TARGET_NAME}.wasm\\n\")
      message(STATUS \"Wrote artifact hash ${TARGET_NAME}.wasm.sha256: \${HS_WASM_SHA256}\")
    ")
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
