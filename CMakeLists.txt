cmake_minimum_required(VERSION 3.10)
project(HolosphereWasm)

# Ensure C++20 for consistency with manual build
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source files
set(SOURCES wasm_bridge.cpp)

# Target name and output name
set(TARGET_NAME holosphere_wasm)

add_executable(${TARGET_NAME} ${SOURCES})

# Include the current directory for all headers
target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Common Emscripten options
target_link_options(${TARGET_NAME} PRIVATE 
    "-lembind"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sEXPORT_NAME=createHolosphereModule"
    "-sSTACK_SIZE=8388608"
)

# Build-type specific flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Configuring Debug build (-O0 -g)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O0" "-g")
    target_link_options(${TARGET_NAME} PRIVATE "-O0" "-g" "-sASSERTIONS=1")
else()
    message(STATUS "Configuring Release build (-O3)")
    target_compile_options(${TARGET_NAME} PRIVATE "-O3")
    target_link_options(${TARGET_NAME} PRIVATE "-O3")
endif()

# Emscripten outputs a .js file as the "executable"
set_target_properties(${TARGET_NAME} PROPERTIES SUFFIX ".js")

# Installation logic
# By default, install to ../daydream relative to the source tree
if(NOT CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # If user provided a prefix, use it
else()
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH "Install location" FORCE)
endif()

# We want to copy both the .js and the generated .wasm file
install(TARGETS ${TARGET_NAME} DESTINATION daydream)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.wasm" DESTINATION daydream)


# --- TESTS ---
add_subdirectory(tests)
