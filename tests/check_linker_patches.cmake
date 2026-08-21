# Pin the routing changes tools/phantasm.ld makes to the stock Teensy 4 linker
# script, and the framework release it was forked from. The script is a verbatim
# copy of framework-arduinoteensy's imxrt1062.ld plus those changes; nothing in
# the build re-derives it, so a framework bump keeps the old stock layout with
# no complaint. Losing a routing change surfaces only as a Teensy size-gate red
# far from its cause -- ITCM past its bank ceiling, or reaction_graph out of
# FLASH -- and only where PlatformIO is installed.
# -D args: SCRIPT (tools/phantasm.ld), PLATFORMIO_INI (platformio.ini).

# Script mode inherits no policies from the project, so every policy would
# otherwise default to OLD. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

file(READ "${SCRIPT}" _script)
file(READ "${PLATFORMIO_INI}" _ini)

set(_missing "")

if(_ini MATCHES "framework-arduinoteensy@([0-9]+\\.[0-9]+\\.[0-9]+)")
  set(_framework "${CMAKE_MATCH_1}")
  if(NOT _script MATCHES "FORKED FROM: framework-arduinoteensy ${_framework}")
    list(APPEND _missing
      "phantasm.ld does not declare FORKED FROM: framework-arduinoteensy ${_framework}, the release platformio.ini pins")
  endif()
else()
  list(APPEND _missing
    "platformio.ini pins no framework-arduinoteensy release to fork from")
endif()

string(FIND "${_script}" ".text.code :" _code_at)
string(FIND "${_script}" ".text.progmem :" _progmem_at)
string(FIND "${_script}" ".text.itcm :" _itcm_at)

if(_code_at EQUAL -1 OR _progmem_at EQUAL -1 OR _itcm_at EQUAL -1)
  list(APPEND _missing
    "phantasm.ld no longer defines all of .text.code, .text.progmem, .text.itcm")
else()
  # Both flash sections must precede .text.itcm's *(.text*) catch-all, which
  # otherwise claims their input sections for ITCM first.
  if(NOT _code_at LESS _itcm_at OR NOT _progmem_at LESS _itcm_at)
    list(APPEND _missing
      ".text.code / .text.progmem no longer precede .text.itcm")
  endif()
  foreach(_section IN ITEMS ".text.code" ".text.progmem")
    if(NOT _script MATCHES "\\${_section} :[^}]*}[ \t]*> FLASH")
      list(APPEND _missing "${_section} is no longer routed to FLASH")
    endif()
  endforeach()

  # Cold library code and the per-effect factories leave ITCM.
  set(_flash_text
    "*libm.a:*(.text .text.*)"
    "*libc_nano.a:*(.text .text.*)"
    "*libstdc++_nano.a:*(.text .text.*)"
    "*libg_nano.a:*(.text .text.*)"
    "*libgcc.a:*(.text .text.*)"
    "*(.text.*construct_effect*)")
  foreach(_spec IN LISTS _flash_text)
    string(FIND "${_script}" "${_spec}" _spec_at)
    if(_spec_at LESS _code_at OR NOT _spec_at LESS _itcm_at)
      list(APPEND _missing "${_spec} is no longer routed by .text.code")
    endif()
  endforeach()

  # Const rodata leaves DTCM.
  set(_rodata "*(SORT_BY_ALIGNMENT(SORT_BY_NAME(.rodata*)))")
  string(FIND "${_script}" "${_rodata}" _rodata_at)
  if(_rodata_at LESS _progmem_at OR NOT _rodata_at LESS _itcm_at)
    list(APPEND _missing "${_rodata} is no longer routed by .text.progmem")
  endif()
endif()

# The stock script's .data claims *(.rodata*) too; leaving that line in place
# would double-claim what .text.progmem now routes to flash.
if(_script MATCHES "\\.data :[^}]*}")
  string(FIND "${CMAKE_MATCH_0}" ".rodata" _rodata_in_data)
  if(NOT _rodata_in_data EQUAL -1)
    list(APPEND _missing ".data claims .rodata back from .text.progmem")
  endif()
else()
  list(APPEND _missing "phantasm.ld no longer defines a .data output section")
endif()

if(_missing)
  string(REPLACE ";" "\n  - " _report "${_missing}")
  message(FATAL_ERROR
    "${SCRIPT} has drifted from the stock script it forks:\n  - "
    "${_report}\n"
    "Re-derive it from framework-arduinoteensy's cores/teensy4/imxrt1062.ld, "
    "re-apply the routing changes its header records, and update the FORKED "
    "FROM line.")
endif()

message(STATUS "phantasm.ld: fork base and both routing changes present")
