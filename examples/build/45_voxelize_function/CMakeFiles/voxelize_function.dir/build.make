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
include 45_voxelize_function/CMakeFiles/voxelize_function.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include 45_voxelize_function/CMakeFiles/voxelize_function.dir/compiler_depend.make

# Include the progress variables for this target.
include 45_voxelize_function/CMakeFiles/voxelize_function.dir/progress.make

# Include the compile flags for this target's objects.
include 45_voxelize_function/CMakeFiles/voxelize_function.dir/flags.make

45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.o: 45_voxelize_function/CMakeFiles/voxelize_function.dir/flags.make
45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.o: /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/45_voxelize_function/main.cpp
45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.o: 45_voxelize_function/CMakeFiles/voxelize_function.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object 45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.o"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/45_voxelize_function && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT 45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.o -MF CMakeFiles/voxelize_function.dir/main.cpp.o.d -o CMakeFiles/voxelize_function.dir/main.cpp.o -c /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/45_voxelize_function/main.cpp

45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/voxelize_function.dir/main.cpp.i"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/45_voxelize_function && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/45_voxelize_function/main.cpp > CMakeFiles/voxelize_function.dir/main.cpp.i

45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/voxelize_function.dir/main.cpp.s"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/45_voxelize_function && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/45_voxelize_function/main.cpp -o CMakeFiles/voxelize_function.dir/main.cpp.s

# Object files for target voxelize_function
voxelize_function_OBJECTS = \
"CMakeFiles/voxelize_function.dir/main.cpp.o"

# External object files for target voxelize_function
voxelize_function_EXTERNAL_OBJECTS =

/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: 45_voxelize_function/CMakeFiles/voxelize_function.dir/main.cpp.o
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: 45_voxelize_function/CMakeFiles/voxelize_function.dir/build.make
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/System/Library/Frameworks/OpenGL.framework
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: imgui/libimgui.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: imgui/glfw/src/libglfw3.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: STB/libSTB.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: _deps/tetgen-build/libtetgen.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: _deps/triangle-build/libtriangle.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: shewchuk_predicates/libshewchuk_predicates.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: _deps/indirect_predicates-build/libindirectPredicates.a
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOImport-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOExport-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOGeometry-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOImage-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOXML-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkRenderingContext2D-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkRenderingFreeType-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkfreetype-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/usr/lib/libz.tbd
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkRenderingVtkJS-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkRenderingSceneGraph-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libjsoncpp.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkImagingCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOLegacy-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkRenderingCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOXMLParser-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkIOCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkFiltersCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkCommonExecutionModel-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkCommonDataModel-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkCommonMisc-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkCommonTransforms-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkCommonMath-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkkissfft-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtkCommonCore-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: /opt/homebrew/lib/libvtksys-9.2.9.2.6.dylib
/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function: 45_voxelize_function/CMakeFiles/voxelize_function.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function"
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/45_voxelize_function && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/voxelize_function.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
45_voxelize_function/CMakeFiles/voxelize_function.dir/build: /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/bin/voxelize_function
.PHONY : 45_voxelize_function/CMakeFiles/voxelize_function.dir/build

45_voxelize_function/CMakeFiles/voxelize_function.dir/clean:
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/45_voxelize_function && $(CMAKE_COMMAND) -P CMakeFiles/voxelize_function.dir/cmake_clean.cmake
.PHONY : 45_voxelize_function/CMakeFiles/voxelize_function.dir/clean

45_voxelize_function/CMakeFiles/voxelize_function.dir/depend:
	cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/45_voxelize_function /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/45_voxelize_function /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/examples/build/45_voxelize_function/CMakeFiles/voxelize_function.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : 45_voxelize_function/CMakeFiles/voxelize_function.dir/depend

