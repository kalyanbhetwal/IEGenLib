# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/apps/cmake/bin/cmake

# The command to remove a file.
RM = /usr/local/apps/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/KALYANBHETWAL/IEGenLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote

# Include any dependencies generated for this target.
include tutorial/CMakeFiles/macro_testing.dir/depend.make

# Include the progress variables for this target.
include tutorial/CMakeFiles/macro_testing.dir/progress.make

# Include the compile flags for this target's objects.
include tutorial/CMakeFiles/macro_testing.dir/flags.make

tutorial/CMakeFiles/macro_testing.dir/macro_testing.cc.o: tutorial/CMakeFiles/macro_testing.dir/flags.make
tutorial/CMakeFiles/macro_testing.dir/macro_testing.cc.o: ../tutorial/macro_testing.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tutorial/CMakeFiles/macro_testing.dir/macro_testing.cc.o"
	cd /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/tutorial && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/macro_testing.dir/macro_testing.cc.o -c /home/KALYANBHETWAL/IEGenLib/tutorial/macro_testing.cc

tutorial/CMakeFiles/macro_testing.dir/macro_testing.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/macro_testing.dir/macro_testing.cc.i"
	cd /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/tutorial && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/KALYANBHETWAL/IEGenLib/tutorial/macro_testing.cc > CMakeFiles/macro_testing.dir/macro_testing.cc.i

tutorial/CMakeFiles/macro_testing.dir/macro_testing.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/macro_testing.dir/macro_testing.cc.s"
	cd /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/tutorial && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/KALYANBHETWAL/IEGenLib/tutorial/macro_testing.cc -o CMakeFiles/macro_testing.dir/macro_testing.cc.s

# Object files for target macro_testing
macro_testing_OBJECTS = \
"CMakeFiles/macro_testing.dir/macro_testing.cc.o"

# External object files for target macro_testing
macro_testing_EXTERNAL_OBJECTS =

bin/tutorial/macro_testing: tutorial/CMakeFiles/macro_testing.dir/macro_testing.cc.o
bin/tutorial/macro_testing: tutorial/CMakeFiles/macro_testing.dir/build.make
bin/tutorial/macro_testing: src/libiegenlib.a
bin/tutorial/macro_testing: tutorial/CMakeFiles/macro_testing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/tutorial/macro_testing"
	cd /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/tutorial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/macro_testing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tutorial/CMakeFiles/macro_testing.dir/build: bin/tutorial/macro_testing

.PHONY : tutorial/CMakeFiles/macro_testing.dir/build

tutorial/CMakeFiles/macro_testing.dir/clean:
	cd /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/tutorial && $(CMAKE_COMMAND) -P CMakeFiles/macro_testing.dir/cmake_clean.cmake
.PHONY : tutorial/CMakeFiles/macro_testing.dir/clean

tutorial/CMakeFiles/macro_testing.dir/depend:
	cd /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/KALYANBHETWAL/IEGenLib /home/KALYANBHETWAL/IEGenLib/tutorial /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/tutorial /home/KALYANBHETWAL/IEGenLib/cmake-build-debug-remote/tutorial/CMakeFiles/macro_testing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tutorial/CMakeFiles/macro_testing.dir/depend

