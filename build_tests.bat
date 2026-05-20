@echo off
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
if %errorlevel% neq 0 (
    echo vcvars64 failed with error %errorlevel%
    exit /b %errorlevel%
)
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
