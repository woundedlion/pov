# Run the relax_bake_gen harness and diff the header tools/relax_bakes.py emits
# from its dump against the committed core/mesh/relax_bakes_generated.h in full.
# unit_relax_bake_verify pins the payload VALUES bit-exact; this pins the file's
# FORM (banner, chunking, declaration layout), so a legitimate regeneration can
# never arrive buried in an emitter reformat.
# -D args: PYTHON_EXE, SCRIPT, HARNESS, DUMP.

# Script mode inherits no policies from the project, so every policy would
# otherwise default to OLD, and the cmake_language(EXIT) below is a 3.29
# feature. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

# Through a file rather than a pipe: execute_process cannot chain into a second
# command without hiding the harness's own exit status.
execute_process(
  COMMAND "${HARNESS}"
  OUTPUT_FILE "${DUMP}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "relax_bake_gen failed (${_rc}):\n${_err}")
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" "${SCRIPT}" check --dump "${DUMP}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "relax_bakes.py check failed (${_rc}):\n${_err}")
endif()

message(STATUS "relax_bake form pin: header text matches the emitter")
