# Run the relax_bake_gen harness and diff the header tools/relax_bakes.py emits
# from its dump against the committed core/mesh/relax_bakes_generated.h in full.
# unit_relax_bake_verify pins the payload VALUES bit-exact; this pins the file's
# FORM (banner, chunking, declaration layout), so a legitimate regeneration can
# never arrive buried in an emitter reformat.
# -D args: PYTHON_EXE, SCRIPT, HARNESS, DUMP.

# Match the top-level CMake policy version in script mode.
cmake_minimum_required(VERSION 3.29)

# Store the dump between the separately checked harness and emitter commands.
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
