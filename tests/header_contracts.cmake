# Compile each test header standalone and twice-included.
file(GLOB_RECURSE _hs_contract_headers CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
set(_hs_contract_stamps)
foreach(_header IN LISTS _hs_contract_headers)
  file(RELATIVE_PATH _relative "${HS_ROOT}" "${_header}")
  string(SHA256 _key "${_relative}")
  set(_source "${CMAKE_CURRENT_BINARY_DIR}/header_contracts/${_key}.cpp")
  set(_stamp "${CMAKE_CURRENT_BINARY_DIR}/header_contracts/${_key}.stamp")
  set(_deps "${CMAKE_CURRENT_BINARY_DIR}/header_contracts/${_key}.d")
  file(CONFIGURE OUTPUT "${_source}" CONTENT "#include \"${_relative}\"\n#include \"${_relative}\"\n")
  add_custom_command(OUTPUT "${_stamp}"
    COMMAND "${CMAKE_CXX_COMPILER}" -std=c++20 -fsyntax-only
      -Wall -Wextra -Werror -ferror-limit=3
      "-I${HS_ROOT}" "-I${CMAKE_CURRENT_BINARY_DIR}"
      "$<$<BOOL:$<TARGET_PROPERTY:run_tests,COMPILE_DEFINITIONS>>:-D$<JOIN:$<TARGET_PROPERTY:run_tests,COMPILE_DEFINITIONS>,;-D>>"
      -MMD -MF "${_deps}" -MT "${_stamp}"
      "${_source}"
    COMMAND "${CMAKE_COMMAND}" -E touch "${_stamp}"
    DEPENDS "${_source}" "${_header}" "${_guard_output}"
    DEPFILE "${_deps}"
    COMMAND_EXPAND_LISTS VERBATIM)
  list(APPEND _hs_contract_stamps "${_stamp}")
endforeach()
add_custom_target(test_header_contracts DEPENDS ${_hs_contract_stamps})
add_dependencies(test_header_contracts death_guard_census)
