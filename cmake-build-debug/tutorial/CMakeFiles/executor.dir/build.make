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
include tutorial/CMakeFiles/executor.dir/depend.make
# Include the progress variables for this target.
include tutorial/CMakeFiles/executor.dir/progress.make

# Include the compile flags for this target's objects.
include tutorial/CMakeFiles/executor.dir/flags.make

tutorial/CMakeFiles/executor.dir/executor.cc.o: tutorial/CMakeFiles/executor.dir/flags.make
tutorial/CMakeFiles/executor.dir/executor.cc.o: ../tutorial/executor.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tutorial/CMakeFiles/executor.dir/executor.cc.o"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/executor.dir/executor.cc.o -c /Users/ant/Documents/AdaptLab/IEGenLib/tutorial/executor.cc

tutorial/CMakeFiles/executor.dir/executor.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/executor.dir/executor.cc.i"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ant/Documents/AdaptLab/IEGenLib/tutorial/executor.cc > CMakeFiles/executor.dir/executor.cc.i

tutorial/CMakeFiles/executor.dir/executor.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/executor.dir/executor.cc.s"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ant/Documents/AdaptLab/IEGenLib/tutorial/executor.cc -o CMakeFiles/executor.dir/executor.cc.s

# Object files for target executor
executor_OBJECTS = \
"CMakeFiles/executor.dir/executor.cc.o"

# External object files for target executor
executor_EXTERNAL_OBJECTS =

bin/tutorial/executor: tutorial/CMakeFiles/executor.dir/executor.cc.o
bin/tutorial/executor: tutorial/CMakeFiles/executor.dir/build.make
bin/tutorial/executor: src/libiegenlib.a
bin/tutorial/executor: tutorial/CMakeFiles/executor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/tutorial/executor"
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/executor.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tutorial/CMakeFiles/executor.dir/build: bin/tutorial/executor
.PHONY : tutorial/CMakeFiles/executor.dir/build

tutorial/CMakeFiles/executor.dir/clean:
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial && $(CMAKE_COMMAND) -P CMakeFiles/executor.dir/cmake_clean.cmake
.PHONY : tutorial/CMakeFiles/executor.dir/clean

tutorial/CMakeFiles/executor.dir/depend:
	cd /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ant/Documents/AdaptLab/IEGenLib /Users/ant/Documents/AdaptLab/IEGenLib/tutorial /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial /Users/ant/Documents/AdaptLab/IEGenLib/cmake-build-debug/tutorial/CMakeFiles/executor.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tutorial/CMakeFiles/executor.dir/depend

