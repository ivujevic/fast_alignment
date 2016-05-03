SET(BOOST_MIN_VERSION "1.60.0")

set(Boost_NO_SYSTEM_PATHS OFF)
set(BOOST_INCLUDEDIR "./boost/include/")
set(BOOST_LIBRARYDIR "./boost/lib/")

find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS serialization system thread program_options REQUIRED)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_library(tac_deps deps.cpp)

target_include_directories(tac_deps SYSTEM INTERFACE
        ${Boost_INCLUDE_DIRS})


target_link_libraries(tac_deps PUBLIC ${Boost_LIBRARIES})





