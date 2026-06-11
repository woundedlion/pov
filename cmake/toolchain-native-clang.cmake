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
    set(CMAKE_C_COMPILER   clang   CACHE FILEPATH "")
    set(CMAKE_CXX_COMPILER clang++ CACHE FILEPATH "")
  endif()
endif()

if(WIN32)
  # CMake's Windows-Clang support links via lld-link, but the emsdk ships only
  # lld.exe. lld is a multicall binary: copied to lld-link.exe it acts as the
  # COFF linker (and embeds manifests itself, so — unlike MSVC link.exe — it
  # does not need rc.exe on PATH at link time). Create the alias next to clang
  # so it is found at both configure and build time.
  if(_hs_clang_dir AND EXISTS "${_hs_clang_dir}/lld.exe"
     AND NOT EXISTS "${_hs_clang_dir}/lld-link.exe")
    file(COPY_FILE "${_hs_clang_dir}/lld.exe" "${_hs_clang_dir}/lld-link.exe")
  endif()
  set(CMAKE_LINKER_TYPE LLD)

  # CMake still enables the RC language for Windows-Clang and test-compiles a
  # stub .rc, so a resource compiler must exist even though we compile no .rc.
  if(NOT CMAKE_RC_COMPILER)
    file(GLOB _hs_rc_candidates
         "$ENV{WindowsSdkVerBinPath}x64/rc.exe"
         "C:/Program Files (x86)/Windows Kits/10/bin/*/x64/rc.exe"
         "C:/Program Files/Windows Kits/10/bin/*/x64/rc.exe")
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
