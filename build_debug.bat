@echo off
set "PATH=%PATH%;C:\Program Files\CMake\bin"
call c:\work\emsdk\emsdk_env.bat
if not exist build_debug mkdir build_debug
cd build_debug
echo Running CMake...
call emcmake cmake -DCMAKE_BUILD_TYPE=Debug -G Ninja ..
if %errorlevel% neq 0 (
    echo CMake failed with error %errorlevel%
    cd ..
    exit /b %errorlevel%
)
echo Running Ninja...
echo Executing: ninja -j 16 -v
ninja -j 16 -v
if %errorlevel% neq 0 (
    echo Ninja failed with error %errorlevel%
    cd ..
    exit /b %errorlevel%
)
echo Installing...
ninja install
if %errorlevel% neq 0 (
    echo Install failed with error %errorlevel%
    cd ..
    exit /b %errorlevel%
)
echo Build completed successfully.
cd ..