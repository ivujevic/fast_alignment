set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/..")

add_library(tac_core
			Tachyon/Base.cpp
			Utils/Bioinformatics/Constants.cpp
			Tachyon/DatabaseElement.cpp
			Utils/Bioinformatics/Evalue.cpp
			Utils/Bioinformatics/Alignment/Opal.cpp
			Utils/Bioinformatics/ScoreMatrix.cpp
			Utils/Bioinformatics/Seg.cpp
			Tachyon/Tachyon.cpp
			Utils/util.cpp Utils/Bioinformatics/Alignment.cpp Utils/Bioinformatics/Alignment.hpp Utils/Bioinformatics/Writer.cpp Utils/Bioinformatics/Writer.hpp Utils/json.hpp
			Utils/Ref.hpp)

target_include_directories(tac_core PUBLIC "{CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(tac_core PRIVATE "{CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(tac_core PUBLIC ${Boost_LIBRARIES})

add_executable(tachyon main.cpp)
target_link_libraries(tachyon tac_core)
