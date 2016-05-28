# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /home/aguo/gitsource/CombBLAS-15_C/CombBLAS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aguo/gitsource/CombBLAS-15_C/CombBLAS

# Include any dependencies generated for this target.
include ReleaseTests/CMakeFiles/MultTest.dir/depend.make

# Include the progress variables for this target.
include ReleaseTests/CMakeFiles/MultTest.dir/progress.make

# Include the compile flags for this target's objects.
include ReleaseTests/CMakeFiles/MultTest.dir/flags.make

ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o: ReleaseTests/CMakeFiles/MultTest.dir/flags.make
ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o: ReleaseTests/MultTest.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && /usr/bin/mpicxx   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/MultTest.dir/MultTest.o -c /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/MultTest.cpp

ReleaseTests/CMakeFiles/MultTest.dir/MultTest.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MultTest.dir/MultTest.i"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/MultTest.cpp > CMakeFiles/MultTest.dir/MultTest.i

ReleaseTests/CMakeFiles/MultTest.dir/MultTest.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MultTest.dir/MultTest.s"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/MultTest.cpp -o CMakeFiles/MultTest.dir/MultTest.s

ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.requires:
.PHONY : ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.requires

ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.provides: ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.requires
	$(MAKE) -f ReleaseTests/CMakeFiles/MultTest.dir/build.make ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.provides.build
.PHONY : ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.provides

ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.provides.build: ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o

# Object files for target MultTest
MultTest_OBJECTS = \
"CMakeFiles/MultTest.dir/MultTest.o"

# External object files for target MultTest
MultTest_EXTERNAL_OBJECTS =

ReleaseTests/MultTest: ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o
ReleaseTests/MultTest: ReleaseTests/CMakeFiles/MultTest.dir/build.make
ReleaseTests/MultTest: libCommGridlib.a
ReleaseTests/MultTest: libMPITypelib.a
ReleaseTests/MultTest: libMemoryPoollib.a
ReleaseTests/MultTest: libHashlib.a
ReleaseTests/MultTest: ReleaseTests/CMakeFiles/MultTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable MultTest"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MultTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ReleaseTests/CMakeFiles/MultTest.dir/build: ReleaseTests/MultTest
.PHONY : ReleaseTests/CMakeFiles/MultTest.dir/build

ReleaseTests/CMakeFiles/MultTest.dir/requires: ReleaseTests/CMakeFiles/MultTest.dir/MultTest.o.requires
.PHONY : ReleaseTests/CMakeFiles/MultTest.dir/requires

ReleaseTests/CMakeFiles/MultTest.dir/clean:
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && $(CMAKE_COMMAND) -P CMakeFiles/MultTest.dir/cmake_clean.cmake
.PHONY : ReleaseTests/CMakeFiles/MultTest.dir/clean

ReleaseTests/CMakeFiles/MultTest.dir/depend:
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aguo/gitsource/CombBLAS-15_C/CombBLAS /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests /home/aguo/gitsource/CombBLAS-15_C/CombBLAS /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/CMakeFiles/MultTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ReleaseTests/CMakeFiles/MultTest.dir/depend

