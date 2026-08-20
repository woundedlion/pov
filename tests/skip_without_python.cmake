# Stands in for a CTest whose producer cannot be built without a Python
# interpreter, so the pair reports as skipped rather than going unregistered: an
# absent test shrinks `ctest`'s total as silently as a passing one raises it.
# -D args: TEST_NAME, SKIP_CODE.

# Script mode inherits no policies from the project, and cmake_language(EXIT) is
# a 3.29 feature. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

message(STATUS "${TEST_NAME}: no python3/python found; skipping")
cmake_language(EXIT ${SKIP_CODE})
