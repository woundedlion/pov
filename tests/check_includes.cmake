# Pin run_tests.cpp's per-module include block against its HS_TEST_MODULE_LIST
# roster. Every roster row references its run_*_tests function, so a missing
# include is already a compile error and the include set is a superset of the
# roster headers. Equal counts therefore prove there is no orphaned include left
# behind after a module was removed. Convention: every tests/test_*.h include in
# run_tests.cpp is a roster module.
# -D args: SRC (path to run_tests.cpp).

file(READ "${SRC}" _text)

string(REGEX MATCHALL "#include \"tests/test_[A-Za-z0-9_]+\\.h\"" _includes "${_text}")
string(REGEX MATCHALL "X\\(\"[A-Za-z0-9_]+\"" _rows "${_text}")

list(LENGTH _includes _ninc)
list(LENGTH _rows _nrow)

if(NOT _ninc EQUAL _nrow)
  message(FATAL_ERROR
    "run_tests.cpp: ${_ninc} test-module includes vs ${_nrow} HS_TEST_MODULE_LIST "
    "rows. An include without a matching roster row compiles dead test source "
    "silently; add the row or drop the orphaned include.")
endif()

message(STATUS "run_tests include pin: ${_ninc} includes match the roster")
