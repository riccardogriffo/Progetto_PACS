# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /u/sw/pkgs/toolchains/gcc-glibc/9/base/bin/cmake

# The command to remove a file.
RM = /u/sw/pkgs/toolchains/gcc-glibc/9/base/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/riccardo/Desktop/codice pacs"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/riccardo/Desktop/codice pacs/build"

# Include any dependencies generated for this target.
include tests/mesh_tests/CMakeFiles/runMeshTest.dir/depend.make

# Include the progress variables for this target.
include tests/mesh_tests/CMakeFiles/runMeshTest.dir/progress.make

# Include the compile flags for this target's objects.
include tests/mesh_tests/CMakeFiles/runMeshTest.dir/flags.make

tests/mesh_tests/CMakeFiles/runMeshTest.dir/MeshTest.cpp.o: tests/mesh_tests/CMakeFiles/runMeshTest.dir/flags.make
tests/mesh_tests/CMakeFiles/runMeshTest.dir/MeshTest.cpp.o: ../tests/mesh_tests/MeshTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/mesh_tests/CMakeFiles/runMeshTest.dir/MeshTest.cpp.o"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runMeshTest.dir/MeshTest.cpp.o -c "/home/riccardo/Desktop/codice pacs/tests/mesh_tests/MeshTest.cpp"

tests/mesh_tests/CMakeFiles/runMeshTest.dir/MeshTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runMeshTest.dir/MeshTest.cpp.i"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/riccardo/Desktop/codice pacs/tests/mesh_tests/MeshTest.cpp" > CMakeFiles/runMeshTest.dir/MeshTest.cpp.i

tests/mesh_tests/CMakeFiles/runMeshTest.dir/MeshTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runMeshTest.dir/MeshTest.cpp.s"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/riccardo/Desktop/codice pacs/tests/mesh_tests/MeshTest.cpp" -o CMakeFiles/runMeshTest.dir/MeshTest.cpp.s

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.o: tests/mesh_tests/CMakeFiles/runMeshTest.dir/flags.make
tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.o: ../src/domain/QuadMesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.o"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.o -c "/home/riccardo/Desktop/codice pacs/src/domain/QuadMesh.cpp"

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.i"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/riccardo/Desktop/codice pacs/src/domain/QuadMesh.cpp" > CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.i

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.s"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/riccardo/Desktop/codice pacs/src/domain/QuadMesh.cpp" -o CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.s

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.o: tests/mesh_tests/CMakeFiles/runMeshTest.dir/flags.make
tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.o: ../src/domain/Node.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.o"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.o -c "/home/riccardo/Desktop/codice pacs/src/domain/Node.cpp"

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.i"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/riccardo/Desktop/codice pacs/src/domain/Node.cpp" > CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.i

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.s"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/riccardo/Desktop/codice pacs/src/domain/Node.cpp" -o CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.s

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.o: tests/mesh_tests/CMakeFiles/runMeshTest.dir/flags.make
tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.o: ../src/domain/Edge.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.o"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.o -c "/home/riccardo/Desktop/codice pacs/src/domain/Edge.cpp"

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.i"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/riccardo/Desktop/codice pacs/src/domain/Edge.cpp" > CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.i

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.s"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/riccardo/Desktop/codice pacs/src/domain/Edge.cpp" -o CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.s

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.o: tests/mesh_tests/CMakeFiles/runMeshTest.dir/flags.make
tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.o: ../src/domain/QuadElement.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.o"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.o -c "/home/riccardo/Desktop/codice pacs/src/domain/QuadElement.cpp"

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.i"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/riccardo/Desktop/codice pacs/src/domain/QuadElement.cpp" > CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.i

tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.s"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/riccardo/Desktop/codice pacs/src/domain/QuadElement.cpp" -o CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.s

# Object files for target runMeshTest
runMeshTest_OBJECTS = \
"CMakeFiles/runMeshTest.dir/MeshTest.cpp.o" \
"CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.o" \
"CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.o" \
"CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.o" \
"CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.o"

# External object files for target runMeshTest
runMeshTest_EXTERNAL_OBJECTS =

tests/mesh_tests/runMeshTest: tests/mesh_tests/CMakeFiles/runMeshTest.dir/MeshTest.cpp.o
tests/mesh_tests/runMeshTest: tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadMesh.cpp.o
tests/mesh_tests/runMeshTest: tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Node.cpp.o
tests/mesh_tests/runMeshTest: tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/Edge.cpp.o
tests/mesh_tests/runMeshTest: tests/mesh_tests/CMakeFiles/runMeshTest.dir/__/__/src/domain/QuadElement.cpp.o
tests/mesh_tests/runMeshTest: tests/mesh_tests/CMakeFiles/runMeshTest.dir/build.make
tests/mesh_tests/runMeshTest: tests/lib/googletest-master/googlemock/gtest/libgtest.a
tests/mesh_tests/runMeshTest: tests/lib/googletest-master/googlemock/gtest/libgtest_main.a
tests/mesh_tests/runMeshTest: tests/lib/googletest-master/googlemock/gtest/libgtest.a
tests/mesh_tests/runMeshTest: tests/mesh_tests/CMakeFiles/runMeshTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable runMeshTest"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/runMeshTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/mesh_tests/CMakeFiles/runMeshTest.dir/build: tests/mesh_tests/runMeshTest

.PHONY : tests/mesh_tests/CMakeFiles/runMeshTest.dir/build

tests/mesh_tests/CMakeFiles/runMeshTest.dir/clean:
	cd "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" && $(CMAKE_COMMAND) -P CMakeFiles/runMeshTest.dir/cmake_clean.cmake
.PHONY : tests/mesh_tests/CMakeFiles/runMeshTest.dir/clean

tests/mesh_tests/CMakeFiles/runMeshTest.dir/depend:
	cd "/home/riccardo/Desktop/codice pacs/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/riccardo/Desktop/codice pacs" "/home/riccardo/Desktop/codice pacs/tests/mesh_tests" "/home/riccardo/Desktop/codice pacs/build" "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests" "/home/riccardo/Desktop/codice pacs/build/tests/mesh_tests/CMakeFiles/runMeshTest.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : tests/mesh_tests/CMakeFiles/runMeshTest.dir/depend

