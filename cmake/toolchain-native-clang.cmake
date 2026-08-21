# Native (non-Emscripten) toolchain for the Holosphere unit tests.
#
# The engine relies on GCC/Clang extensions (__attribute__((always_inline)),
# noinline, etc.) that MSVC does not accept, so the native test build must use
# Clang. On Windows this also smooths over two CMake/Clang papercuts so the
# suite configures from a plain shell (no Visual Studio Developer Prompt):
#   * provides the lld-link the emsdk omits, and
#   * locates the Windows SDK resource compiler.
#
# Compiler resolution order:
#   1. An explicitly provided CMAKE_CXX_COMPILER (respected, not overridden).
#   2. $ENV{EMSDK}/upstream/bin/clang++         (set by emsdk_env).
#   3. <repo>/../emsdk/upstream/bin/clang++     (sibling emsdk checkout).
#   4. clang/clang++ from PATH                  (Linux/macOS CI, or PATH setups).

# --- Locate the Clang bin directory (used for both the compiler and lld) ---
set(_hs_clang_dir "")
if(DEFINED ENV{EMSDK} AND EXISTS "$ENV{EMSDK}/upstream/bin")
  set(_hs_clang_dir "$ENV{EMSDK}/upstream/bin")
else()
  # cmake/ lives at the repo root; ../.. is the parent dir (sibling emsdk).
  get_filename_component(_hs_repo_parent "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
  if(EXISTS "${_hs_repo_parent}/emsdk/upstream/bin")
    set(_hs_clang_dir "${_hs_repo_parent}/emsdk/upstream/bin")
  endif()
endif()

if(NOT CMAKE_CXX_COMPILER)
  if(_hs_clang_dir AND WIN32)
    set(CMAKE_C_COMPILER   "${_hs_clang_dir}/clang.exe"   CACHE FILEPATH "")
    set(CMAKE_CXX_COMPILER "${_hs_clang_dir}/clang++.exe" CACHE FILEPATH "")
  elseif(_hs_clang_dir)
    set(CMAKE_C_COMPILER   "${_hs_clang_dir}/clang"   CACHE FILEPATH "")
    set(CMAKE_CXX_COMPILER "${_hs_clang_dir}/clang++" CACHE FILEPATH "")
  else()
    # Whatever clang PATH resolves; tests/CMakeLists.txt holds the version floor.
    set(CMAKE_C_COMPILER   clang   CACHE FILEPATH "")
    set(CMAKE_CXX_COMPILER clang++ CACHE FILEPATH "")
  endif()
endif()

if(WIN32)
  # CMake's Windows-Clang support links via lld-link, but the emsdk ships only
  # lld.exe. lld is a multicall binary: invoked as lld-link.exe it acts as the
  # COFF linker (and embeds manifests itself, so — unlike MSVC link.exe — it
  # does not need rc.exe on PATH at link time). Stage the alias in the BUILD
  # TREE (always writable) rather than writing it into the emsdk install, which
  # hard-fails on a read-only toolchain checkout; add the build tree to clang's
  # program-search prefix (-B) so it finds lld-link there at link time without
  # the alias sitting next to clang or on PATH. (The MSVC-target driver ignores
  # --ld-path for lld-link, so -B is the mechanism that actually applies.)
  if(_hs_clang_dir AND EXISTS "${_hs_clang_dir}/lld.exe")
    # Unconditional copy: an existing build dir must not keep a stale alias from
    # a previous emsdk when the clang it links through has been upgraded.
    file(COPY_FILE "${_hs_clang_dir}/lld.exe" "${CMAKE_BINARY_DIR}/lld-link.exe")
    string(APPEND CMAKE_EXE_LINKER_FLAGS_INIT " -B\"${CMAKE_BINARY_DIR}\"")
  endif()
  set(CMAKE_LINKER_TYPE LLD)

  # CMake still enables the RC language for Windows-Clang and test-compiles a
  # stub .rc, so a resource compiler must exist even though we compile no .rc.
  if(NOT CMAKE_RC_COMPILER)
    # Windows Kits installs under a Program Files root, so take those roots from
    # the environment; the literal C: paths are the last resort for a stripped
    # environment.
    set(_hs_sdk_roots "")
    foreach(_hs_pf "$ENV{ProgramFiles\(x86\)}" "$ENV{ProgramFiles}")
      if(NOT _hs_pf STREQUAL "")
        file(TO_CMAKE_PATH "${_hs_pf}" _hs_pf_slashed)
        list(APPEND _hs_sdk_roots "${_hs_pf_slashed}")
      endif()
    endforeach()
    if(NOT "$ENV{SystemDrive}" STREQUAL "")
      list(APPEND _hs_sdk_roots "$ENV{SystemDrive}/Program Files (x86)"
                                "$ENV{SystemDrive}/Program Files")
    endif()
    list(APPEND _hs_sdk_roots "C:/Program Files (x86)" "C:/Program Files")
    list(REMOVE_DUPLICATES _hs_sdk_roots)

    set(_hs_rc_globs "$ENV{WindowsSdkVerBinPath}x64/rc.exe")
    foreach(_hs_sdk_root IN LISTS _hs_sdk_roots)
      list(APPEND _hs_rc_globs "${_hs_sdk_root}/Windows Kits/10/bin/*/x64/rc.exe")
    endforeach()
    file(GLOB _hs_rc_candidates ${_hs_rc_globs})
    if(_hs_rc_candidates)
      # NATURAL so version components sort numerically (10.0.22621 > 10.0.9...);
      # a plain lexical sort would rank 10.0.9xxxx above 10.0.22xxx and pick an
      # older SDK's rc.exe.
      list(SORT _hs_rc_candidates COMPARE NATURAL)
      list(GET _hs_rc_candidates -1 _hs_rc)  # highest SDK version sorts last
      set(CMAKE_RC_COMPILER "${_hs_rc}" CACHE FILEPATH "")
    endif()
  endif()
endif()
