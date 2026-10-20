# Pin tools/phantasm.ld routing changes and the source Teensy framework release.
# -D args: SCRIPT (tools/phantasm.ld), PLATFORMIO_INI (platformio.ini).

# Match the top-level CMake policy version in script mode.
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

string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _script "${_script}")

string(FIND "${_script}" ".text.code :" _code_at)
string(FIND "${_script}" ".text.progmem :" _progmem_at)
string(FIND "${_script}" ".text.itcm :" _itcm_at)

if(_code_at EQUAL -1 OR _progmem_at EQUAL -1 OR _itcm_at EQUAL -1)
  list(APPEND _missing
    "phantasm.ld no longer defines all of .text.code, .text.progmem, .text.itcm")
else()
  # .text.code, then .text.progmem, both before .text.itcm's *(.text*) catch-all.
  if(NOT _code_at LESS _progmem_at OR NOT _progmem_at LESS _itcm_at)
    list(APPEND _missing
      ".text.code, .text.progmem, .text.itcm are no longer in that order")
  endif()
  foreach(_section IN ITEMS ".text.code" ".text.progmem")
    if(NOT _script MATCHES "\\${_section} :[^}]*}[ \t]*> FLASH")
      list(APPEND _missing "${_section} is no longer routed to FLASH")
    endif()
  endforeach()

  foreach(_hot_route IN ITEMS
      "*(.text._ZN4Plot18walk_rigid_samplesI*)"
      "*(.text._ZN4Plot15walk_dense_runsI*)")
    string(FIND "${_script}" "${_hot_route}" _hot_at)
    if(_hot_at LESS _itcm_at)
      list(APPEND _missing "timed hot route missing from ITCM: ${_hot_route}")
    endif()
  endforeach()

  # Cold library code and the per-effect factories leave ITCM.
  set(_flash_text
    "*libm.a:*(.text .text.*)"
    "*libc_nano.a:*(.text .text.*)"
    "*libstdc++_nano.a:*(.text .text.*)"
    "*libg_nano.a:*(.text .text.*)"
    "*libgcc.a:*(.text .text.*)"
    "*libFrameworkArduino.a:usb.c.o(.text .text.*)"
    "*libFrameworkArduino.a:usb_serial.c.o(.text .text.*)"
    "*(.text.*construct_effect*)"
    "*(.text.*construct_profiled*)"
    "*(.text._ZN8SPIClass16beginTransactionE*)"
    "*(.text._ZN2hs12CycleCounter8log_nodeE*)"
    "*(.text._ZN4math6detail15upsample_framesE*)"
    "*(.text._ZZN12POVSegmentedI*8run_showI*)"
    "*(.text.*ProfiledEffect*10draw_frameEv*)"
    "*(.text._ZN6Effect18update_render_clipEv*)"
    "*(.text._ZN6Effect8set_clipE*)"
    "*(.text._ZN9Animation6MotionI*4stepER6Canvas)"
    "*(.text._ZNK9Animation6MotionI*10path_frameEf)"
    "*(.text._ZN9Animation8RotationI*4stepER6Canvas)"
    "*(.text._ZN9Animation10RandomWalkI*4stepER6Canvas)"
    "*(.text._ZN9Animation26PausableParamAnimationBaseI*4stepER6Canvas)"
    "*(.text._ZN9Animation14MobiusWarpBaseI*4stepER6Canvas)"
    "*(.text._ZN9Animation18MobiusWarpEvolving4stepER6Canvas)"
    "*(.text._ZN9Animation6Ripple4stepER6Canvas)"
    "*(.text._ZN9Animation5Noise4stepER6Canvas)"
    "*(.text._ZN9Animation12NoiseProduct4stepER6Canvas)"
    "*(.text._ZN9Animation8BallDropI*4stepER6Canvas)"
    "*(.text._ZN3SDF4FaceC1E*)"
    "*(.text._ZN3SDF4FaceC2E*)"
    "*(.text._ZN3SDF4Face19compute_full_boundsE*)"
    "*(.text._ZN5Arena4makeIN8Pullback14ComposedEffectI*5StateE*)"
    "*(.text._ZN8Pullback14ComposedEffectI*5StateC1Ev)"
    "*(.text._ZN8Pullback14ComposedEffectI*5StateC2Ev)"
    "*(.text._ZN8Pullback15TemporalHistoryI*C1Ev)"
    "*(.text._ZN8Pullback15TemporalHistoryI*C2Ev)")
  foreach(_spec IN LISTS _flash_text)
    string(FIND "${_script}" "${_spec}" _spec_at)
    if(_spec_at LESS _code_at OR NOT _spec_at LESS _progmem_at)
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

# .data must not also claim the .rodata that .text.progmem routes to flash.
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
