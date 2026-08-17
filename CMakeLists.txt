# 3.29 is the real floor: CMakePresets schema v6 needs >=3.25 (3.21-3.24 reject
# the file outright) and the native test toolchain's CMAKE_LINKER_TYPE was added
# in 3.29 (silently ignored on 3.25-3.28, so the Windows tests build falls back
# to the wrong linker — the exact failure toolchain-native-clang.cmake exists to
# prevent). Declare it honestly so an old CMake fails loudly instead of mis-linking.
cmake_minimum_required(VERSION 3.29)
project(Holosphere CXX)

set(HS_PULLBACK_CAPTURE_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH
    "Checkout whose sources a pullback capture build compiles")
set(HS_PULLBACK_CAPTURE_PRESET "" CACHE STRING
    "CMake preset attested by pullback capture provenance")

# Two build modes share this file, selected by the toolchain:
#   * Emscripten toolchain  -> WASM simulator module (daydream)   [EMSCRIPTEN]
#   * Native clang toolchain -> unit-test executable + CTest      [else]
# Use the presets in CMakePresets.json (cmake --preset wasm-release / tests).

# Enable ccache if available (both modes benefit). find_program stores the
# resolved executable PATH (or <var>-NOTFOUND), so use CMake's CCACHE_PROGRAM
# convention rather than a _FOUND name that reads as a boolean but holds a path.
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
  message(STATUS "ccache enabled: ${CCACHE_PROGRAM}")
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

  # Shipping module: sources come from CMAKE_CURRENT_SOURCE_DIR, the root the
  # install rules below record provenance against. Only the capture backends
  # follow HS_PULLBACK_CAPTURE_SOURCE_ROOT.
  set(TARGET_NAME holosphere_wasm)
  add_executable(${TARGET_NAME}
    ${CMAKE_CURRENT_SOURCE_DIR}/targets/wasm/wasm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/core/engine/memory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/core/engine/static_storage.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/core/engine/reaction_graph.cpp)

  target_include_directories(${TARGET_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core)

  # Warning hygiene for the JS bridge: targets/wasm/ is built by no other target,
  # so without this it is the one first-party TU with zero warning coverage. Not
  # -Werror: unlike the native suite, third-party emscripten system headers ride
  # this compile line. The gate is in CI instead — the wasm job greps both build
  # logs for warnings under first-party paths and fails on any
  # (.github/workflows/ci.yml), the same first-party/third-party split the device
  # firmware ratchet uses (tools/teensy_warnings.py).
  target_compile_options(${TARGET_NAME} PRIVATE "-Wall" "-Wextra")

  target_link_options(${TARGET_NAME} PRIVATE
    "-lembind"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sEXPORT_NAME=createHolosphereModule"
    # No runtime code generation, so the pages that load this module need only
    # 'wasm-unsafe-eval' in their CSP rather than 'unsafe-eval'. embind normally
    # crafts one invoker per binding with `new Function` during module creation;
    # DYNAMIC_EXECUTION=0 removes that, and EMBIND_AOT emits the same invokers
    # into the glue at link time instead. Both are needed: DYNAMIC_EXECUTION=0
    # alone falls back to a generic reflective invoker roughly 8x slower per
    # bound call. Cost of EMBIND_AOT is ~5.5 KB of glue and ~2 s of link time.
    "-sDYNAMIC_EXECUTION=0"
    "-sEMBIND_AOT=1")
  # -sSTACK_SIZE is deliberately NOT set here: each build type sets it exactly
  # once in its own block below (8 KB release / 64 KB debug). Setting it in this
  # shared block too would make the effective value depend on link-line ordering
  # between two -sSTACK_SIZE flags (last one wins), so a future reorder of these
  # blocks could silently revert debug to 8 KB.

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WASM Debug build (-O0 -g -sASSERTIONS, 64 KB stack)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O0" "-g")
    # 64 KB stack for debug only. -O0 disables inlining and stack-slot
    # coalescing, inflating frame sizes severalfold over the -O3 release build,
    # so the tight release budget overflows under no optimization. CI does
    # runtime-smoke this debug build, so size it generously to hold the inflated
    # -O0 frames rather than rely on -sASSERTIONS to merely report the abort.
    target_link_options(${TARGET_NAME} PRIVATE "-O0" "-g" "-sASSERTIONS=1"
      "-sSTACK_SIZE=65536")
  elseif(CMAKE_BUILD_TYPE STREQUAL "PullbackStrictFP")
    message(STATUS "WASM pullback strict-FP diagnostic build (-O3 -fno-fast-math -ffp-contract=off, no LTO)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O3" "-fno-fast-math"
      "-ffp-contract=off" "-msimd128")
    target_link_options(${TARGET_NAME} PRIVATE "-O3" "-fno-fast-math"
      "-ffp-contract=off" "-msimd128" "-sSTACK_SIZE=8192")
  else()
    # -fno-finite-math-only is required AFTER -ffast-math: -ffast-math implies
    # -ffinite-math-only, which folds std::isfinite() to constant true and lets
    # the compiler assume no NaN/Inf module-wide. The render sink relies on real
    # finite semantics (filter.h: a NaN coord must not silently UB-cast and
    # bypass the clip/wrap guards). We keep the rest of -ffast-math (reassoc,
    # FMA, reciprocal, -fno-math-errno) for the SIMD float throughput.
    # The device firmware build carries the SAME -ffast-math -fno-finite-math-only
    # pair (platformio.ini, shared [env].build_flags) so both targets share float
    # reassoc/FMA semantics. Note this does NOT make float results bit-identical
    # across ARM and WASM (different ISA/libm/FMA contraction) — the sim/device
    # parity contract is the integer wrap guards, not floats.
    message(STATUS "WASM Release build (-O3 -ffast-math -fno-finite-math-only -msimd128 -flto, 8 KB stack)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O3" "-ffast-math" "-fno-finite-math-only" "-msimd128" "-flto")
    # The -ffast-math / -fno-finite-math-only pair also rides the link line. With
    # -flto the math semantics live in per-function IR attributes that the LTO
    # link consumes (it does not re-derive them from the command line), so the
    # link-line copy is a no-op for today's toolchain — but it closes the gap if a
    # future link-time math pass ever keys off the link flags instead.
    # 8 KB stack: the tight production budget, sized against the inlined -O3 build.
    target_link_options(${TARGET_NAME} PRIVATE "-O3" "-ffast-math"
      "-fno-finite-math-only" "-msimd128" "-flto" "-sSTACK_SIZE=8192")
  endif()

  # Optional dev-only WASM bindings (off by default, so a normal build is
  # unchanged). Enables MeshOps::getMaxBounds() — a roster measurement tool for
  # sizing MAX_VERTS-style constants, with no UI/runtime consumer (see
  # HS_WASM_DEV_BINDINGS in targets/wasm/mesh_ops_bindings.h). Enable with:
  #   cmake --preset wasm-release -DHS_WASM_DEV_BINDINGS=ON
  option(HS_WASM_DEV_BINDINGS "Re-export dev-only WASM bindings (getMaxBounds)" OFF)
  if(HS_WASM_DEV_BINDINGS)
    target_compile_definitions(${TARGET_NAME} PRIVATE HS_WASM_DEV_BINDINGS)
    message(STATUS "WASM dev bindings enabled (HS_WASM_DEV_BINDINGS): getMaxBounds re-exported")
    set(HS_WASM_DEV_BINDINGS_STATE ON)
  else()
    set(HS_WASM_DEV_BINDINGS_STATE OFF)
  endif()

  # Emscripten emits a .js loader as the "executable".
  set_target_properties(${TARGET_NAME} PROPERTIES SUFFIX ".js")

  # Out of the default build: its own LTO link costs as much as the module's and
  # nothing but a capture run consumes the .cjs, so tools/pullback_crosscheck.py
  # names it as an explicit --target. Its four sources are compiled by the module
  # target and by tests/'s pullback_capture_native, so excluding it hides no
  # source from CI.
  add_executable(pullback_capture_wasm EXCLUDE_FROM_ALL
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/pullback_capture_native.cpp
    ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/engine/memory.cpp
    ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/engine/static_storage.cpp
    ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/engine/reaction_graph.cpp)
  target_include_directories(pullback_capture_wasm PRIVATE
    ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}
    ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core)
  target_compile_definitions(pullback_capture_wasm PRIVATE
    HS_TIMELINE_MAX_ANIM_BYTES=256
    HS_GLOBAL_ARENA_BYTES=8388608
    HS_ENABLE_EFFECT_REGISTRY=1
    HS_ENABLE_PARAM_GUI_BRIDGE=1
    HS_ENABLE_TEST_HOOKS=1
    HS_ENABLE_TEST_ORACLES=1
    HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND=0
    HS_ENABLE_EFFECT_CONTROL_API=1)
  target_compile_options(pullback_capture_wasm PRIVATE "-Wall" "-Wextra")
  target_link_options(pullback_capture_wasm PRIVATE
    "-sALLOW_MEMORY_GROWTH=1"
    "-sENVIRONMENT=node"
    "-sNODERAWFS=1")
  if(CMAKE_BUILD_TYPE STREQUAL "PullbackStrictFP")
    target_compile_options(pullback_capture_wasm PRIVATE "-O3"
      "-fno-fast-math" "-ffp-contract=off" "-msimd128")
    target_link_options(pullback_capture_wasm PRIVATE "-O3"
      "-fno-fast-math" "-ffp-contract=off" "-msimd128"
      "-sSTACK_SIZE=65536")
  else()
    target_compile_options(pullback_capture_wasm PRIVATE "-O3"
      "-ffast-math" "-fno-finite-math-only" "-msimd128" "-flto")
    target_link_options(pullback_capture_wasm PRIVATE "-O3"
      "-ffast-math" "-fno-finite-math-only" "-msimd128" "-flto"
      "-sSTACK_SIZE=65536")
  endif()
  set_target_properties(pullback_capture_wasm PROPERTIES SUFFIX ".cjs")

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

  # The install overwrites README.md and mirrors docs/screenshots by wiping the
  # destination first, so a typo'd or stale cache value pointing at an unrelated
  # directory deletes and overwrites paths there. Require a marker only a
  # daydream checkout carries. daydream.js is tracked there and is written by no
  # install rule here, so a previous run cannot have planted it.
  set(HS_DAYDREAM_MARKER "${DAYDREAM_DIR}/daydream.js")

  if(EXISTS "${HS_DAYDREAM_MARKER}")
    message(STATUS "Simulator module installs to DAYDREAM_DIR: ${DAYDREAM_DIR}")
    install(TARGETS ${TARGET_NAME} DESTINATION "${DAYDREAM_DIR}")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.wasm"
            DESTINATION "${DAYDREAM_DIR}")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/README.md"
            DESTINATION "${DAYDREAM_DIR}")
    # The segment-map golden, emitted from hardware/pov_segment_map.h and pinned
    # against it by the unit_pov_segment_map_golden CTest. The simulator's
    # tests/segment_crosscheck.test.js reads it as the firmware reference, so it
    # travels the same install path as the module rather than being hand-ported
    # into JS. daydream's deploy gate diffs its committed copy against this file in
    # the engine checkout it pins, so a stale copy there fails too.
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/hardware/pov_segment_map.json"
            DESTINATION "${DAYDREAM_DIR}")
    install(FILES
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/shader_workbench.mjs"
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/sha256.mjs"
            DESTINATION "${DAYDREAM_DIR}/shader")
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/patterns/"
            DESTINATION "${DAYDREAM_DIR}/shader/patterns"
            FILES_MATCHING
            PATTERN "*.shader.json"
            PATTERN "shaderball_migration.json")
    # The gallery is mirrored, not merged: the destination is wiped first so a
    # screenshot dropped from the engine stops being served. install(DIRECTORY)
    # installs nothing for an absent or empty source, which would leave the wipe
    # standing as a silent deletion of the live gallery — so the source is
    # checked at install time, immediately before the wipe.
    install(CODE "
      file(GLOB HS_SCREENSHOTS \"${CMAKE_CURRENT_SOURCE_DIR}/docs/screenshots/*.png\")
      if(NOT HS_SCREENSHOTS)
        message(FATAL_ERROR \"No PNGs under ${CMAKE_CURRENT_SOURCE_DIR}/docs/screenshots; refusing to wipe ${DAYDREAM_DIR}/docs/screenshots\")
      endif()
      file(REMOVE_RECURSE \"${DAYDREAM_DIR}/docs/screenshots\")")
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/docs/screenshots/"
            DESTINATION "${DAYDREAM_DIR}/docs/screenshots")

    # Record the engine commit this WASM module was built from, the content hashes
    # of the module's two shipped halves (.wasm binary, .js glue), and the toolchain
    # that produced them, alongside it. The daydream deploy gate reads the SHA to pin
    # its engine checkout to the exact source the committed module corresponds to
    # (instead of a decoupled engine HEAD), and asserts the committed files match the
    # recorded hashes so a stale or swapped binary — or stale glue — can't ride a
    # green gate. All are captured at install time (immediately after the build) from
    # the same build event, and committed together with the artifacts they describe.
    install(CODE "
      execute_process(
        COMMAND git -C \"${CMAKE_CURRENT_SOURCE_DIR}\" rev-parse HEAD
        OUTPUT_VARIABLE HS_ENGINE_SHA
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE HS_GIT_RC
        ERROR_QUIET)
      if(HS_GIT_RC EQUAL 0)
        # Dirty-tree guard: a module built from uncommitted (tracked) source records
        # provenance the deploy gate cannot reproduce — the gate checks out the bare
        # SHA (clean) and tests that, while the shipped binary came from dirty source.
        # Mark the SHA -dirty (git describe convention) so the gate's checkout of that
        # ref fails hard instead of deploying an unverifiable artifact. Untracked files
        # don't count (build output, node_modules); only tracked modifications do.
        # A git that cannot answer counts as dirty: empty output from a failed call
        # is indistinguishable from a clean tree, and recording a bare SHA for it
        # is the unverifiable provenance this block exists to prevent.
        execute_process(
          COMMAND git -C \"${CMAKE_CURRENT_SOURCE_DIR}\" status --porcelain --untracked-files=no
          OUTPUT_VARIABLE HS_GIT_DIRTY
          OUTPUT_STRIP_TRAILING_WHITESPACE
          RESULT_VARIABLE HS_GIT_DIRTY_RC)
        if(NOT HS_GIT_DIRTY_RC EQUAL 0 OR NOT \"\${HS_GIT_DIRTY}\" STREQUAL \"\")
          set(HS_ENGINE_SHA \"\${HS_ENGINE_SHA}-dirty\")
          message(WARNING \"Engine tree has uncommitted changes (or git could not report on it); recording ${TARGET_NAME}.sha as \${HS_ENGINE_SHA}. The deploy gate will reject this ref — rebuild from a clean (committed) tree before deploying.\")
        endif()
        file(CONFIGURE OUTPUT \"${DAYDREAM_DIR}/${TARGET_NAME}.sha\"
             CONTENT \"\${HS_ENGINE_SHA}\\n\" @ONLY NEWLINE_STYLE LF)
        message(STATUS \"Wrote engine provenance ${TARGET_NAME}.sha: \${HS_ENGINE_SHA}\")

        # sha256sum -c manifest, one '<hash>  <filename>' line per artifact, checked
        # from the daydream root. Covers both halves of the module: the binary and the
        # Emscripten glue that instantiates it — stale or hand-edited glue (embind
        # signature drift, wrong memory config) white-screens the live page against an
        # otherwise byte-verified binary. Kept under the .wasm.sha256 name the deploy
        # gate already passes to `sha256sum -c`.
        # Written only alongside the .sha so the provenance set is never half-written
        # (a .wasm.sha256 without a .sha is a partial set the deploy gate rejects).
        # Git stores the generated glue and provenance files as LF on every host,
        # and emcc emits CRLF glue on Windows. file(CONFIGURE ... NEWLINE_STYLE LF)
        # is the only writer that lands LF on every host — file(WRITE)/file(APPEND)
        # open in text mode and re-expand LF to CRLF there — but it also expands
        # @VAR@ tokens, and the hash below is taken over its output. Escaping every
        # literal @ through a variable that expands back to @ makes that expansion
        # an identity, so glue bytes reach the hash exactly as emitted.
        set(HS_AT \"@\")
        file(READ \"${DAYDREAM_DIR}/${TARGET_NAME}.js\" HS_GLUE_CONTENT)
        string(REPLACE \"@\" \"@HS_AT@\" HS_GLUE_CONTENT \"\${HS_GLUE_CONTENT}\")
        file(CONFIGURE OUTPUT \"${DAYDREAM_DIR}/${TARGET_NAME}.js\"
             CONTENT \"\${HS_GLUE_CONTENT}\" @ONLY NEWLINE_STYLE LF)
        file(SHA256 \"${DAYDREAM_DIR}/${TARGET_NAME}.wasm\" HS_WASM_SHA256)
        file(SHA256 \"${DAYDREAM_DIR}/${TARGET_NAME}.js\" HS_GLUE_SHA256)
        file(CONFIGURE OUTPUT \"${DAYDREAM_DIR}/${TARGET_NAME}.wasm.sha256\"
             CONTENT \"\${HS_WASM_SHA256}  ${TARGET_NAME}.wasm\\n\${HS_GLUE_SHA256}  ${TARGET_NAME}.js\\n\"
             @ONLY NEWLINE_STYLE LF)
        message(STATUS \"Wrote artifact hashes ${TARGET_NAME}.wasm.sha256: .wasm=\${HS_WASM_SHA256} .js=\${HS_GLUE_SHA256}\")

        # Toolchain identity of the build event. The .sha pins the SOURCE, but the
        # deploy gate re-verifies that source with its own compiler, so a module
        # built by an unpinned local emsdk otherwise records provenance the gate
        # cannot distinguish from CI's pinned one (ci.yml wasm job: emsdk 5.0.0).
        # A separate additive file, not extra lines in the two above: the gate feeds
        # .wasm.sha256 to `sha256sum -c` and requires .sha to be a bare 40-hex ref.
        # Written with the pair above, keeping the provenance set all-or-nothing.
        #
        # build_type and dev_bindings record the CONFIGURATION as well as the
        # compiler: the same source and the same emsdk produce a different module
        # under wasm-debug (-O0, 64 KB stack, assertions) or with
        # -DHS_WASM_DEV_BINDINGS=ON (extra exports). Without them a module built
        # from either rides a gate that only re-verifies source and emsdk.
        file(CONFIGURE OUTPUT \"${DAYDREAM_DIR}/${TARGET_NAME}.toolchain\"
             CONTENT \"emsdk ${EMSCRIPTEN_VERSION}\\nclang ${CMAKE_CXX_COMPILER_VERSION}\\nbuild_type ${CMAKE_BUILD_TYPE}\\ndev_bindings ${HS_WASM_DEV_BINDINGS_STATE}\\n\"
             @ONLY NEWLINE_STYLE LF)
        message(STATUS \"Wrote toolchain provenance ${TARGET_NAME}.toolchain: emsdk ${EMSCRIPTEN_VERSION}, clang ${CMAKE_CXX_COMPILER_VERSION}, ${CMAKE_BUILD_TYPE}, dev bindings ${HS_WASM_DEV_BINDINGS_STATE}\")
      else()
        # Best-effort provenance: git is absent or this is not a working tree (an
        # exported tarball / non-git checkout). Warn and omit ALL markers so a plain
        # 'just install' still succeeds outside a git checkout, rather than the old
        # FATAL_ERROR. The deploy gate then finds no provenance and refuses to deploy
        # an unverifiable artifact — which is correct; deploy only from a git checkout.
        message(WARNING \"git unavailable or not a working tree (rc=\${HS_GIT_RC}); skipping ${TARGET_NAME} provenance markers (.sha / .wasm.sha256 / .toolchain). Install completes, but the deploy gate has nothing to verify — deploy only from a git checkout.\")
      endif()
    ")
  elseif(EXISTS "${DAYDREAM_DIR}")
    # An existing destination without the marker is a wrong path, not a missing
    # checkout: fail the install rather than write into it.
    message(WARNING
      "DAYDREAM_DIR '${DAYDREAM_DIR}' exists but holds no daydream.js — "
      "refusing to install into it.")
    install(CODE "message(FATAL_ERROR \"DAYDREAM_DIR '${DAYDREAM_DIR}' exists but is not a daydream checkout (no daydream.js); refusing to overwrite its README.md and wipe its docs/screenshots. Point -DDAYDREAM_DIR at the daydream checkout and re-run cmake.\")")
  else()
    message(WARNING
      "DAYDREAM_DIR '${DAYDREAM_DIR}' not found — skipping simulator-module "
      "install (set -DDAYDREAM_DIR=<path> and re-run cmake to enable).")
    # The configure-time warning above is not reprinted by a build against an
    # already-configured tree, which is every `just install` after the first.
    install(CODE "message(WARNING \"DAYDREAM_DIR '${DAYDREAM_DIR}' not found; installed nothing. Set -DDAYDREAM_DIR=<path> and re-run cmake.\")")
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
