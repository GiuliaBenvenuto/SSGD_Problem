# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build

# Include any dependencies generated for this target.
include _deps/tetgen-build/CMakeFiles/tetgen.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include _deps/tetgen-build/CMakeFiles/tetgen.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/tetgen-build/CMakeFiles/tetgen.dir/progress.make

# Include the compile flags for this target's objects.
include _deps/tetgen-build/CMakeFiles/tetgen.dir/flags.make

_deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.o: _deps/tetgen-build/CMakeFiles/tetgen.dir/flags.make
_deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.o: _deps/tetgen-src/tetgen.cxx
_deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.o: _deps/tetgen-build/CMakeFiles/tetgen.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object _deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.o"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT _deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.o -MF CMakeFiles/tetgen.dir/tetgen.cxx.o.d -o CMakeFiles/tetgen.dir/tetgen.cxx.o -c /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-src/tetgen.cxx

_deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tetgen.dir/tetgen.cxx.i"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-src/tetgen.cxx > CMakeFiles/tetgen.dir/tetgen.cxx.i

_deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tetgen.dir/tetgen.cxx.s"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-src/tetgen.cxx -o CMakeFiles/tetgen.dir/tetgen.cxx.s

_deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.o: _deps/tetgen-build/CMakeFiles/tetgen.dir/flags.make
_deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.o: _deps/tetgen-src/predicates.cxx
_deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.o: _deps/tetgen-build/CMakeFiles/tetgen.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object _deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.o"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT _deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.o -MF CMakeFiles/tetgen.dir/predicates.cxx.o.d -o CMakeFiles/tetgen.dir/predicates.cxx.o -c /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-src/predicates.cxx

_deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tetgen.dir/predicates.cxx.i"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-src/predicates.cxx > CMakeFiles/tetgen.dir/predicates.cxx.i

_deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tetgen.dir/predicates.cxx.s"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-src/predicates.cxx -o CMakeFiles/tetgen.dir/predicates.cxx.s

# Object files for target tetgen
tetgen_OBJECTS = \
"CMakeFiles/tetgen.dir/tetgen.cxx.o" \
"CMakeFiles/tetgen.dir/predicates.cxx.o"

# External object files for target tetgen
tetgen_EXTERNAL_OBJECTS =

_deps/tetgen-build/libtetgen.a: _deps/tetgen-build/CMakeFiles/tetgen.dir/tetgen.cxx.o
_deps/tetgen-build/libtetgen.a: _deps/tetgen-build/CMakeFiles/tetgen.dir/predicates.cxx.o
_deps/tetgen-build/libtetgen.a: _deps/tetgen-build/CMakeFiles/tetgen.dir/build.make
_deps/tetgen-build/libtetgen.a: _deps/tetgen-build/CMakeFiles/tetgen.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libtetgen.a"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && $(CMAKE_COMMAND) -P CMakeFiles/tetgen.dir/cmake_clean_target.cmake
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tetgen.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
_deps/tetgen-build/CMakeFiles/tetgen.dir/build: _deps/tetgen-build/libtetgen.a
.PHONY : _deps/tetgen-build/CMakeFiles/tetgen.dir/build

_deps/tetgen-build/CMakeFiles/tetgen.dir/clean:
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build && $(CMAKE_COMMAND) -P CMakeFiles/tetgen.dir/cmake_clean.cmake
.PHONY : _deps/tetgen-build/CMakeFiles/tetgen.dir/clean

_deps/tetgen-build/CMakeFiles/tetgen.dir/depend:
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-src /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/_deps/tetgen-build/CMakeFiles/tetgen.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : _deps/tetgen-build/CMakeFiles/tetgen.dir/depend

