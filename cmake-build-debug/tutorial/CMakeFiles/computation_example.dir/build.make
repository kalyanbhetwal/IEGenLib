# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ant/Documents/AdaptLab/IEGenLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug

# Include any dependencies generated for this target.
include tutorial/CMakeFiles/computation_example.dir/depend.make
# Include the progress variables for this target.
include tutorial/CMakeFiles/computation_example.dir/progress.make

# Include the compile flags for this target's objects.
include tutorial/CMakeFiles/computation_example.dir/flags.make

tutorial/CMakeFiles/computation_example.dir/computation_example.cc.o: tutorial/CMakeFiles/computation_example.dir/flags.make
tutorial/CMakeFiles/computation_example.dir/computation_example.cc.o: ../tutorial/computation_example.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tutorial/CMakeFiles/computation_example.dir/computation_example.cc.o"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/computation_example.dir/computation_example.cc.o -c /Users/ant/Documents/AdaptLab/IEGenLib/tutorial/computation_example.cc

tutorial/CMakeFiles/computation_example.dir/computation_example.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/computation_example.dir/computation_example.cc.i"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ant/Documents/AdaptLab/IEGenLib/tutorial/computation_example.cc > CMakeFiles/computation_example.dir/computation_example.cc.i

tutorial/CMakeFiles/computation_example.dir/computation_example.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/computation_example.dir/computation_example.cc.s"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ant/Documents/AdaptLab/IEGenLib/tutorial/computation_example.cc -o CMakeFiles/computation_example.dir/computation_example.cc.s

# Object files for target computation_example
computation_example_OBJECTS = \
"CMakeFiles/computation_example.dir/computation_example.cc.o"

# External object files for target computation_example
computation_example_EXTERNAL_OBJECTS =

bin/tutorial/computation_example: tutorial/CMakeFiles/computation_example.dir/computation_example.cc.o
bin/tutorial/computation_example: tutorial/CMakeFiles/computation_example.dir/build.make
bin/tutorial/computation_example: src/libiegenlib.a
bin/tutorial/computation_example: tutorial/CMakeFiles/computation_example.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/tutorial/computation_example"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/computation_example.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tutorial/CMakeFiles/computation_example.dir/build: bin/tutorial/computation_example
.PHONY : tutorial/CMakeFiles/computation_example.dir/build

tutorial/CMakeFiles/computation_example.dir/clean:
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && $(CMAKE_COMMAND) -P CMakeFiles/computation_example.dir/cmake_clean.cmake
.PHONY : tutorial/CMakeFiles/computation_example.dir/clean

tutorial/CMakeFiles/computation_example.dir/depend:
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ant/Documents/AdaptLab/IEGenLib /Users/ant/Documents/AdaptLab/IEGenLib/tutorial /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial/CMakeFiles/computation_example.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tutorial/CMakeFiles/computation_example.dir/depend

