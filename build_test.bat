
@echo off
set "PATH=%PATH%;C:\Program Files\CMake\bin"
call c:\work\emsdk\emsdk_env.bat
if not exist build_test mkdir build_test
cd build_test
echo Running CMake...
call emcmake cmake -DCMAKE_BUILD_TYPE=Debug -G Ninja ..
if %errorlevel% neq 0 (
    echo CMake failed with error %errorlevel%
    cd ..
    exit /b %errorlevel%
)
echo Running Ninja for test_winding...
ninja test_winding
if %errorlevel% neq 0 (
    echo Ninja failed with error %errorlevel%
    cd ..
    exit /b %errorlevel%
)
echo Build completed successfully.
cd ..
