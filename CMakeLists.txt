cmake_minimum_required(VERSION 3.0)
project(tachyon)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")


include(CodeCoverage)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -g -O3 -std=c++11 -march=native -lpthread -MMD -MP -MF -Wall -Wextra \
						-Wno-sign-compare -Wno-return-type -fno-strict-aliasing")


endif()

enable_testing()
add_subdirectory(deps)
add_subdirectory(src)
