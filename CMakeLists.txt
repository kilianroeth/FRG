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

find_package(HDF5 REQUIRED COMPONENTS C)

target_include_directories(FRG PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(FRG PRIVATE ${HDF5_LIBRARIES})

# including directories
target_include_directories(${PROJECT_NAME} PRIVATE include)

if (MSVC) 
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Release>:-O2>
    )
endif()
