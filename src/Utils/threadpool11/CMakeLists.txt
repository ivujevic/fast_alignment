add_definitions(-Dthreadpool11_EXPORTING)

include_directories(include)

add_library(threadpool11 STATIC src/pool.cpp src/worker.cpp include/threadpool11/pool.hpp include/threadpool11/worker.hpp include/threadpool11/utility.hpp)

if (CMAKE_COMPILER_IS_GNUCXX)
    target_link_libraries(threadpool11 pthread)
endif()

