# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vujevic/TacLast/Tachyon

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vujevic/TacLast/Tachyon/build

# Include any dependencies generated for this target.
include src/CMakeFiles/runquery.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/runquery.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/runquery.dir/flags.make

src/CMakeFiles/runquery.dir/RunQuery.cpp.o: src/CMakeFiles/runquery.dir/flags.make
src/CMakeFiles/runquery.dir/RunQuery.cpp.o: ../src/RunQuery.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vujevic/TacLast/Tachyon/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/runquery.dir/RunQuery.cpp.o"
	cd /home/vujevic/TacLast/Tachyon/build/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/runquery.dir/RunQuery.cpp.o -c /home/vujevic/TacLast/Tachyon/src/RunQuery.cpp

src/CMakeFiles/runquery.dir/RunQuery.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runquery.dir/RunQuery.cpp.i"
	cd /home/vujevic/TacLast/Tachyon/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/vujevic/TacLast/Tachyon/src/RunQuery.cpp > CMakeFiles/runquery.dir/RunQuery.cpp.i

src/CMakeFiles/runquery.dir/RunQuery.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runquery.dir/RunQuery.cpp.s"
	cd /home/vujevic/TacLast/Tachyon/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/vujevic/TacLast/Tachyon/src/RunQuery.cpp -o CMakeFiles/runquery.dir/RunQuery.cpp.s

src/CMakeFiles/runquery.dir/RunQuery.cpp.o.requires:
.PHONY : src/CMakeFiles/runquery.dir/RunQuery.cpp.o.requires

src/CMakeFiles/runquery.dir/RunQuery.cpp.o.provides: src/CMakeFiles/runquery.dir/RunQuery.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/runquery.dir/build.make src/CMakeFiles/runquery.dir/RunQuery.cpp.o.provides.build
.PHONY : src/CMakeFiles/runquery.dir/RunQuery.cpp.o.provides

src/CMakeFiles/runquery.dir/RunQuery.cpp.o.provides.build: src/CMakeFiles/runquery.dir/RunQuery.cpp.o

# Object files for target runquery
runquery_OBJECTS = \
"CMakeFiles/runquery.dir/RunQuery.cpp.o"

# External object files for target runquery
runquery_EXTERNAL_OBJECTS =

runquery: src/CMakeFiles/runquery.dir/RunQuery.cpp.o
runquery: src/CMakeFiles/runquery.dir/build.make
runquery: deps/librq_deps.a
runquery: /usr/lib/x86_64-linux-gnu/libboost_serialization.so
runquery: /usr/lib/x86_64-linux-gnu/libboost_system.so
runquery: /usr/lib/x86_64-linux-gnu/libboost_thread.so
runquery: /usr/lib/x86_64-linux-gnu/libboost_program_options.so
runquery: /usr/lib/x86_64-linux-gnu/libpthread.so
runquery: src/CMakeFiles/runquery.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../runquery"
	cd /home/vujevic/TacLast/Tachyon/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/runquery.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/runquery.dir/build: runquery
.PHONY : src/CMakeFiles/runquery.dir/build

src/CMakeFiles/runquery.dir/requires: src/CMakeFiles/runquery.dir/RunQuery.cpp.o.requires
.PHONY : src/CMakeFiles/runquery.dir/requires

src/CMakeFiles/runquery.dir/clean:
	cd /home/vujevic/TacLast/Tachyon/build/src && $(CMAKE_COMMAND) -P CMakeFiles/runquery.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/runquery.dir/clean

src/CMakeFiles/runquery.dir/depend:
	cd /home/vujevic/TacLast/Tachyon/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vujevic/TacLast/Tachyon /home/vujevic/TacLast/Tachyon/src /home/vujevic/TacLast/Tachyon/build /home/vujevic/TacLast/Tachyon/build/src /home/vujevic/TacLast/Tachyon/build/src/CMakeFiles/runquery.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/runquery.dir/depend
