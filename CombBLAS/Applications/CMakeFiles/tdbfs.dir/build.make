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
CMAKE_SOURCE_DIR = /home/allenzou/git/CombBLAS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/allenzou/git/CombBLAS

# Include any dependencies generated for this target.
include Applications/CMakeFiles/tdbfs.dir/depend.make

# Include the progress variables for this target.
include Applications/CMakeFiles/tdbfs.dir/progress.make

# Include the compile flags for this target's objects.
include Applications/CMakeFiles/tdbfs.dir/flags.make

Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o: Applications/CMakeFiles/tdbfs.dir/flags.make
Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o: Applications/TopDownBFS.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/allenzou/git/CombBLAS/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o"
	cd /home/allenzou/git/CombBLAS/Applications && /usr/bin/mpicxx   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/tdbfs.dir/TopDownBFS.o -c /home/allenzou/git/CombBLAS/Applications/TopDownBFS.cpp

Applications/CMakeFiles/tdbfs.dir/TopDownBFS.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tdbfs.dir/TopDownBFS.i"
	cd /home/allenzou/git/CombBLAS/Applications && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/allenzou/git/CombBLAS/Applications/TopDownBFS.cpp > CMakeFiles/tdbfs.dir/TopDownBFS.i

Applications/CMakeFiles/tdbfs.dir/TopDownBFS.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tdbfs.dir/TopDownBFS.s"
	cd /home/allenzou/git/CombBLAS/Applications && /usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/allenzou/git/CombBLAS/Applications/TopDownBFS.cpp -o CMakeFiles/tdbfs.dir/TopDownBFS.s

Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.requires:
.PHONY : Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.requires

Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.provides: Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.requires
	$(MAKE) -f Applications/CMakeFiles/tdbfs.dir/build.make Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.provides.build
.PHONY : Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.provides

Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.provides.build: Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o

# Object files for target tdbfs
tdbfs_OBJECTS = \
"CMakeFiles/tdbfs.dir/TopDownBFS.o"

# External object files for target tdbfs
tdbfs_EXTERNAL_OBJECTS =

Applications/tdbfs: Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o
Applications/tdbfs: Applications/CMakeFiles/tdbfs.dir/build.make
Applications/tdbfs: libCommGridlib.a
Applications/tdbfs: libMPITypelib.a
Applications/tdbfs: libMemoryPoollib.a
Applications/tdbfs: graph500-1.2/generator/libGraphGenlib.a
Applications/tdbfs: libHashlib.a
Applications/tdbfs: Applications/CMakeFiles/tdbfs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable tdbfs"
	cd /home/allenzou/git/CombBLAS/Applications && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tdbfs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Applications/CMakeFiles/tdbfs.dir/build: Applications/tdbfs
.PHONY : Applications/CMakeFiles/tdbfs.dir/build

Applications/CMakeFiles/tdbfs.dir/requires: Applications/CMakeFiles/tdbfs.dir/TopDownBFS.o.requires
.PHONY : Applications/CMakeFiles/tdbfs.dir/requires

Applications/CMakeFiles/tdbfs.dir/clean:
	cd /home/allenzou/git/CombBLAS/Applications && $(CMAKE_COMMAND) -P CMakeFiles/tdbfs.dir/cmake_clean.cmake
.PHONY : Applications/CMakeFiles/tdbfs.dir/clean

Applications/CMakeFiles/tdbfs.dir/depend:
	cd /home/allenzou/git/CombBLAS && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/allenzou/git/CombBLAS /home/allenzou/git/CombBLAS/Applications /home/allenzou/git/CombBLAS /home/allenzou/git/CombBLAS/Applications /home/allenzou/git/CombBLAS/Applications/CMakeFiles/tdbfs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Applications/CMakeFiles/tdbfs.dir/depend
