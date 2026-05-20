@echo off
if not exist build_test mkdir build_test
echo Compiling tests...
C:\work\emsdk\upstream\bin\clang++.exe -std=c++20 -O0 -g -I . -I core -o build_test\run_tests.exe tests\run_tests.cpp
if %errorlevel% neq 0 (
    echo Compile failed with error %errorlevel%
    exit /b %errorlevel%
)
echo Running tests...
build_test\run_tests.exe
exit /b %errorlevel%
