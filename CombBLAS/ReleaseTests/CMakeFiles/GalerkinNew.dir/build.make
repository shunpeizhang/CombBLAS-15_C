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
include ReleaseTests/CMakeFiles/GalerkinNew.dir/depend.make

# Include the progress variables for this target.
include ReleaseTests/CMakeFiles/GalerkinNew.dir/progress.make

# Include the compile flags for this target's objects.
include ReleaseTests/CMakeFiles/GalerkinNew.dir/flags.make

ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o: ReleaseTests/CMakeFiles/GalerkinNew.dir/flags.make
ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o: ReleaseTests/GalerkinNew.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && /usr/bin/mpicxx   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/GalerkinNew.dir/GalerkinNew.o -c /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/GalerkinNew.cpp

ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GalerkinNew.dir/GalerkinNew.i"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/GalerkinNew.cpp > CMakeFiles/GalerkinNew.dir/GalerkinNew.i

ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GalerkinNew.dir/GalerkinNew.s"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/GalerkinNew.cpp -o CMakeFiles/GalerkinNew.dir/GalerkinNew.s

ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.requires:
.PHONY : ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.requires

ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.provides: ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.requires
	$(MAKE) -f ReleaseTests/CMakeFiles/GalerkinNew.dir/build.make ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.provides.build
.PHONY : ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.provides

ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.provides.build: ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o

# Object files for target GalerkinNew
GalerkinNew_OBJECTS = \
"CMakeFiles/GalerkinNew.dir/GalerkinNew.o"

# External object files for target GalerkinNew
GalerkinNew_EXTERNAL_OBJECTS =

ReleaseTests/GalerkinNew: ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o
ReleaseTests/GalerkinNew: ReleaseTests/CMakeFiles/GalerkinNew.dir/build.make
ReleaseTests/GalerkinNew: libCommGridlib.a
ReleaseTests/GalerkinNew: libMPITypelib.a
ReleaseTests/GalerkinNew: libMemoryPoollib.a
ReleaseTests/GalerkinNew: libHashlib.a
ReleaseTests/GalerkinNew: ReleaseTests/CMakeFiles/GalerkinNew.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable GalerkinNew"
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GalerkinNew.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ReleaseTests/CMakeFiles/GalerkinNew.dir/build: ReleaseTests/GalerkinNew
.PHONY : ReleaseTests/CMakeFiles/GalerkinNew.dir/build

ReleaseTests/CMakeFiles/GalerkinNew.dir/requires: ReleaseTests/CMakeFiles/GalerkinNew.dir/GalerkinNew.o.requires
.PHONY : ReleaseTests/CMakeFiles/GalerkinNew.dir/requires

ReleaseTests/CMakeFiles/GalerkinNew.dir/clean:
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests && $(CMAKE_COMMAND) -P CMakeFiles/GalerkinNew.dir/cmake_clean.cmake
.PHONY : ReleaseTests/CMakeFiles/GalerkinNew.dir/clean

ReleaseTests/CMakeFiles/GalerkinNew.dir/depend:
	cd /home/aguo/gitsource/CombBLAS-15_C/CombBLAS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aguo/gitsource/CombBLAS-15_C/CombBLAS /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests /home/aguo/gitsource/CombBLAS-15_C/CombBLAS /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests /home/aguo/gitsource/CombBLAS-15_C/CombBLAS/ReleaseTests/CMakeFiles/GalerkinNew.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ReleaseTests/CMakeFiles/GalerkinNew.dir/depend

