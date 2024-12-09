cmake_minimum_required(VERSION 3.22)

project(fluid
    VERSION 1.0
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(TYPES "" CACHE STRING "Types for compilation")

message(STATUS ${TYPES})

add_compile_definitions(TYPES=${TYPES})

add_compile_options("-Wall")
add_compile_options("-Werror")
add_compile_options("-Wextra")
add_compile_options("-DDEBUG")

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE REPO_SRC src/*.cpp)

add_executable(fluid src/fluid.cpp)
