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
include 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/compiler_depend.make

# Include the progress variables for this target.
include 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/progress.make

# Include the compile flags for this target's objects.
include 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/flags.make

13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o: 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/flags.make
13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o: /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/13_polyhedral_mesh_generation/main.cpp
13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o: 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/13_polyhedral_mesh_generation && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o -MF CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o.d -o CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o -c /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/13_polyhedral_mesh_generation/main.cpp

13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.i"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/13_polyhedral_mesh_generation && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/13_polyhedral_mesh_generation/main.cpp > CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.i

13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.s"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/13_polyhedral_mesh_generation && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/13_polyhedral_mesh_generation/main.cpp -o CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.s

# Object files for target polyhedral_mesh_generation
polyhedral_mesh_generation_OBJECTS = \
"CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o"

# External object files for target polyhedral_mesh_generation
polyhedral_mesh_generation_EXTERNAL_OBJECTS =

/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/main.cpp.o
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/build.make
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/System/Library/Frameworks/OpenGL.framework
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: imgui/libimgui.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: imgui/glfw/src/libglfw3.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: STB/libSTB.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: _deps/tetgen-build/libtetgen.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: _deps/triangle-build/libtriangle.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: shewchuk_predicates/libshewchuk_predicates.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: _deps/indirect_predicates-build/libindirectPredicates.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOImport-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOExport-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOGeometry-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOImage-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOXML-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkRenderingContext2D-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkRenderingFreeType-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkfreetype-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/usr/lib/libz.tbd
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkRenderingVtkJS-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkRenderingSceneGraph-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libjsoncpp.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkImagingCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOLegacy-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkRenderingCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOXMLParser-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkIOCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkFiltersCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkCommonExecutionModel-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkCommonDataModel-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkCommonMisc-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkCommonTransforms-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkCommonMath-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkkissfft-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtkCommonCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: /opt/homebrew/lib/libvtksys-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation: 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/13_polyhedral_mesh_generation && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/polyhedral_mesh_generation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/build: /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/polyhedral_mesh_generation
.PHONY : 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/build

13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/clean:
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/13_polyhedral_mesh_generation && $(CMAKE_COMMAND) -P CMakeFiles/polyhedral_mesh_generation.dir/cmake_clean.cmake
.PHONY : 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/clean

13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/depend:
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/13_polyhedral_mesh_generation /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/13_polyhedral_mesh_generation /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : 13_polyhedral_mesh_generation/CMakeFiles/polyhedral_mesh_generation.dir/depend

