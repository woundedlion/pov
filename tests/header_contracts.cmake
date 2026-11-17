# Compile each test header standalone and twice-included.
file(GLOB_RECURSE _hs_contract_headers CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
add_custom_target(test_header_contracts)
foreach(_header IN LISTS _hs_contract_headers)
  file(RELATIVE_PATH _relative "${HS_ROOT}" "${_header}")
  string(SHA256 _key "${_relative}")
  set(_source "${CMAKE_CURRENT_BINARY_DIR}/header_contracts/${_key}.cpp")
  file(CONFIGURE OUTPUT "${_source}" CONTENT "#include \"${_relative}\"\n#include \"${_relative}\"\n")
  add_custom_target(hs_header_${_key}
    COMMAND "${CMAKE_CXX_COMPILER}" -std=c++20 -fsyntax-only
      -Wall -Wextra -Werror -ferror-limit=3
      "-I${HS_ROOT}" "-I${CMAKE_CURRENT_BINARY_DIR}"
      "$<$<BOOL:$<TARGET_PROPERTY:run_tests,COMPILE_DEFINITIONS>>:-D$<JOIN:$<TARGET_PROPERTY:run_tests,COMPILE_DEFINITIONS>,;-D>>"
      "${_source}"
    COMMAND_EXPAND_LISTS VERBATIM)
  add_dependencies(hs_header_${_key} death_guard_census)
  add_dependencies(test_header_contracts hs_header_${_key})
endforeach()
