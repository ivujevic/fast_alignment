#set(Boost_NO_SYSTEM_PATHS ON)
#set(Boost_USE_STATIC_LIBS OFF)
#set(BOOST_INCLUDEDIR "")
#set(BOOST_LIBRARYDIR "")

find_package(Boost COMPONENTS serialization system thread program_options REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_library(tac_deps deps.cpp)
target_include_directories(tac_deps SYSTEM INTERFACE
        ${Boost_INCLUDE_DIRS})


target_link_libraries(tac_deps PUBLIC ${Boost_LIBRARIES})
