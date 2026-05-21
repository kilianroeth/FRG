cmake_minimum_required(VERSION 3.10)

# set the project name
project(FRG VERSION 1.0 LANGUAGES CXX)
# export compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)

# collect all source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# add exectuable
add_executable(${PROJECT_NAME} ${SOURCES})

# including directories
target_include_directories(${PROJECT_NAME} PRIVATE include)

target_compile_options(${PROJECT_NAME} PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Release>:-O2>
)
