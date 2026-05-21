@echo off
set "PATH=%PATH%;C:\Program Files\CMake\bin"
call c:\work\emsdk\emsdk_env.bat
if not exist build_release mkdir build_release
cd build_release
echo Running CMake Release...
call emcmake cmake -DCMAKE_BUILD_TYPE=Release -G Ninja ..
if %errorlevel% neq 0 (
    cd ..
    echo CMake failed
    exit /b %errorlevel%
)
echo Running Ninja...
echo Executing: ninja -j 16 -v
ninja -j 16 -v
if %errorlevel% neq 0 (
    cd ..
    echo Ninja failed
    exit /b %errorlevel%
)
echo Installing...
ninja install
if %errorlevel% neq 0 (
    cd ..
    echo Install failed
    exit /b %errorlevel%
)
cd ..
echo Release Build completed successfully.
