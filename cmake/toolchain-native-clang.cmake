# Native (non-Emscripten) Clang toolchain for the Holosphere unit tests; the
# engine uses GCC/Clang extensions MSVC rejects. On Windows it also provides
# lld-link and the SDK resource compiler, so no Developer Prompt is needed.
#
# Compiler resolution order: an explicit CMAKE_CXX_COMPILER, then
# $ENV{EMSDK}/upstream/bin, then a sibling <repo>/../emsdk, then PATH.

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
    # Whatever clang PATH resolves.
    set(CMAKE_C_COMPILER   clang   CACHE FILEPATH "")
    set(CMAKE_CXX_COMPILER clang++ CACHE FILEPATH "")
  endif()
endif()

if(WIN32)
  # Stage emsdk's lld.exe as the COFF multicall alias lld-link.exe; -B searches
  # the build tree for it.
  if(_hs_clang_dir AND EXISTS "${_hs_clang_dir}/lld.exe")
    file(COPY_FILE "${_hs_clang_dir}/lld.exe" "${CMAKE_BINARY_DIR}/lld-link.exe"
         ONLY_IF_DIFFERENT)
    string(APPEND CMAKE_EXE_LINKER_FLAGS_INIT " -B\"${CMAKE_BINARY_DIR}\"")
  endif()
  set(CMAKE_LINKER_TYPE LLD)

  # CMake test-compiles a stub .rc for Windows-Clang, so a resource compiler
  # must exist.
  if(NOT CMAKE_RC_COMPILER)
    # Program Files roots from the environment; literal C: paths as fallback.
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
      # NATURAL sorts version components numerically (10.0.22621 > 10.0.9...).
      list(SORT _hs_rc_candidates COMPARE NATURAL)
      list(GET _hs_rc_candidates -1 _hs_rc)  # highest SDK version sorts last
      set(CMAKE_RC_COMPILER "${_hs_rc}" CACHE FILEPATH "")
    endif()
  endif()
endif()
