set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/..")
include(XcodeSettings)

add_library(tac_core
		Base.cpp
		Constants.cpp
		DatabaseElement.cpp
		Evalue.cpp
		opal.cpp
		ScoreMatrix.cpp
		Seg.cpp
		tachyon.cpp
		util.cpp alignment.cpp alignment.hpp writer.cpp writer.hpp)

target_include_directories(tac_core PUBLIC "{CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(tac_core PRIVATE "{CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(tac_core PUBLIC tac_deps)
adjust_xcode_settings(tac_core)

add_executable(tachyon main.cpp)
target_link_libraries(tachyon tac_core)
adjust_xcode_settings(tachyon)


add_executable(runquery runquery.cpp)
target_link_libraries(runquery tac_deps)

