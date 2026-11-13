cmake_minimum_required(VERSION 3.29)

file(READ "${SCRIPT}" baseline)
file(MAKE_DIRECTORY "${TEST_ROOT}")

function(check_routing name source diagnostic)
  set(path "${TEST_ROOT}/${name}.ld")
  file(WRITE "${path}" "${source}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" "-DSCRIPT=${path}"
      "-DPLATFORMIO_INI=${PLATFORMIO_INI}"
      -P "${CMAKE_CURRENT_LIST_DIR}/check_linker_patches.cmake"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(diagnostic STREQUAL "")
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "clean routing rejected: ${output}${error}")
    endif()
  elseif(result EQUAL 0 OR NOT "${output}${error}" MATCHES "${diagnostic}")
    message(FATAL_ERROR "${name} did not refuse the expected routing: ${output}${error}")
  endif()
endfunction()

check_routing(baseline "${baseline}" "")
string(REPLACE ".data :" ".data : { *(.rodata*) }\n.data :" changed "${baseline}")
check_routing(data_rodata "${changed}" ".data claims .rodata")
string(REPLACE ".text.code :" ".text.swap :" changed "${baseline}")
string(REPLACE ".text.progmem :" ".text.code :" changed "${changed}")
string(REPLACE ".text.swap :" ".text.progmem :" changed "${changed}")
check_routing(section_order "${changed}" "are no longer in that order")
string(REPLACE "INCLUDE phantasm_hot_layout.ld" "" changed "${baseline}")
string(REPLACE "*libm.a:*(.text .text.*)"
  "*libm.a:*(.text .text.*)\nINCLUDE phantasm_hot_layout.ld" changed "${changed}")
check_routing(late_include "${changed}" "hot-code layout INCLUDE must sit")
string(REPLACE "hs_hot_layout_end = .;" "" changed "${baseline}")
string(REPLACE "*libm.a:*(.text .text.*)"
  "*libm.a:*(.text .text.*)\nhs_hot_layout_end = .;" changed "${changed}")
check_routing(late_layout_end "${changed}" "hot-code layout INCLUDE must sit")
string(REPLACE "*libgcc.a:*(.text .text.*)" "" changed "${baseline}")
check_routing(missing_library "${changed}" "libgcc.*no longer routed by .text.code")
string(REPLACE "*(SORT_BY_ALIGNMENT(SORT_BY_NAME(.rodata*)))" "" changed "${baseline}")
check_routing(missing_rodata "${changed}" "no longer routed by .text.progmem")
set(setup_route_index 0)
foreach(route IN ITEMS
    "*(.text._ZN13PaletteCyclerC[12]Ev)"
    "*(.text._ZN8Pullback16OuterCameraStateILb1EEC[12]Ev)"
    "*(.text._ZN5Arena10allocate_nI8ParamDefEEPT_j)")
  string(REPLACE "${route}" "" changed "${baseline}")
  check_routing("missing_setup_${setup_route_index}" "${changed}"
    "no longer routed by .text.code")
  math(EXPR setup_route_index "${setup_route_index} + 1")
endforeach()
message(STATUS "linker routing controls passed")
