cmake_minimum_required(VERSION 3.0)
project(tachyon)

set(CMAKE_VERBOSE_MAKEFILE )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if( NOT CMAKE_BUILD_TYPE)
    message("Setting build type to 'Release' as none one was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++17 -march=native -fPIC -Werror")

if( CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -g -O1")
elseif( CMAKE_BUILD_TYPE MATCHES Release)
    message("Release mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
find_package(Boost COMPONENTS serialization program_options REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
else()
    message("Can not find Boost, please install it.")
    include_directories()
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

add_subdirectory(src)

