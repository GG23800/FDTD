# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /home/jerome/Bureau/Github/FDTD

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jerome/Bureau/Github/FDTD

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jerome/Bureau/Github/FDTD/CMakeFiles /home/jerome/Bureau/Github/FDTD/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jerome/Bureau/Github/FDTD/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named HeatFDTD

# Build rule for target.
HeatFDTD: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 HeatFDTD
.PHONY : HeatFDTD

# fast build rule for target.
HeatFDTD/fast:
	$(MAKE) -f CMakeFiles/HeatFDTD.dir/build.make CMakeFiles/HeatFDTD.dir/build
.PHONY : HeatFDTD/fast

#=============================================================================
# Target rules for targets named fdtd_test

# Build rule for target.
fdtd_test: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 fdtd_test
.PHONY : fdtd_test

# fast build rule for target.
fdtd_test/fast:
	$(MAKE) -f test/CMakeFiles/fdtd_test.dir/build.make test/CMakeFiles/fdtd_test.dir/build
.PHONY : fdtd_test/fast

src/heat_fdtd_1d.o: src/heat_fdtd_1d.cpp.o

.PHONY : src/heat_fdtd_1d.o

# target to build an object file
src/heat_fdtd_1d.cpp.o:
	$(MAKE) -f CMakeFiles/HeatFDTD.dir/build.make CMakeFiles/HeatFDTD.dir/src/heat_fdtd_1d.cpp.o
.PHONY : src/heat_fdtd_1d.cpp.o

src/heat_fdtd_1d.i: src/heat_fdtd_1d.cpp.i

.PHONY : src/heat_fdtd_1d.i

# target to preprocess a source file
src/heat_fdtd_1d.cpp.i:
	$(MAKE) -f CMakeFiles/HeatFDTD.dir/build.make CMakeFiles/HeatFDTD.dir/src/heat_fdtd_1d.cpp.i
.PHONY : src/heat_fdtd_1d.cpp.i

src/heat_fdtd_1d.s: src/heat_fdtd_1d.cpp.s

.PHONY : src/heat_fdtd_1d.s

# target to generate assembly for a file
src/heat_fdtd_1d.cpp.s:
	$(MAKE) -f CMakeFiles/HeatFDTD.dir/build.make CMakeFiles/HeatFDTD.dir/src/heat_fdtd_1d.cpp.s
.PHONY : src/heat_fdtd_1d.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... HeatFDTD"
	@echo "... edit_cache"
	@echo "... fdtd_test"
	@echo "... src/heat_fdtd_1d.o"
	@echo "... src/heat_fdtd_1d.i"
	@echo "... src/heat_fdtd_1d.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

