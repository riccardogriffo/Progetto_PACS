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
include tests/example_tests/CMakeFiles/runSimpleTest.dir/depend.make

# Include the progress variables for this target.
include tests/example_tests/CMakeFiles/runSimpleTest.dir/progress.make

# Include the compile flags for this target's objects.
include tests/example_tests/CMakeFiles/runSimpleTest.dir/flags.make

tests/example_tests/CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.o: tests/example_tests/CMakeFiles/runSimpleTest.dir/flags.make
tests/example_tests/CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.o: ../tests/example_tests/SimpleTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/example_tests/CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.o"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/example_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.o -c "/home/riccardo/Desktop/codice pacs/tests/example_tests/SimpleTest.cpp"

tests/example_tests/CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.i"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/example_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/riccardo/Desktop/codice pacs/tests/example_tests/SimpleTest.cpp" > CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.i

tests/example_tests/CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.s"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/example_tests" && /u/sw/pkgs/toolchains/gcc-glibc/9/prefix/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/riccardo/Desktop/codice pacs/tests/example_tests/SimpleTest.cpp" -o CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.s

# Object files for target runSimpleTest
runSimpleTest_OBJECTS = \
"CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.o"

# External object files for target runSimpleTest
runSimpleTest_EXTERNAL_OBJECTS =

tests/example_tests/runSimpleTest: tests/example_tests/CMakeFiles/runSimpleTest.dir/SimpleTest.cpp.o
tests/example_tests/runSimpleTest: tests/example_tests/CMakeFiles/runSimpleTest.dir/build.make
tests/example_tests/runSimpleTest: tests/lib/googletest-master/googlemock/gtest/libgtest.a
tests/example_tests/runSimpleTest: tests/lib/googletest-master/googlemock/gtest/libgtest_main.a
tests/example_tests/runSimpleTest: tests/lib/googletest-master/googlemock/gtest/libgtest.a
tests/example_tests/runSimpleTest: tests/example_tests/CMakeFiles/runSimpleTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/riccardo/Desktop/codice pacs/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable runSimpleTest"
	cd "/home/riccardo/Desktop/codice pacs/build/tests/example_tests" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/runSimpleTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/example_tests/CMakeFiles/runSimpleTest.dir/build: tests/example_tests/runSimpleTest

.PHONY : tests/example_tests/CMakeFiles/runSimpleTest.dir/build

tests/example_tests/CMakeFiles/runSimpleTest.dir/clean:
	cd "/home/riccardo/Desktop/codice pacs/build/tests/example_tests" && $(CMAKE_COMMAND) -P CMakeFiles/runSimpleTest.dir/cmake_clean.cmake
.PHONY : tests/example_tests/CMakeFiles/runSimpleTest.dir/clean

tests/example_tests/CMakeFiles/runSimpleTest.dir/depend:
	cd "/home/riccardo/Desktop/codice pacs/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/riccardo/Desktop/codice pacs" "/home/riccardo/Desktop/codice pacs/tests/example_tests" "/home/riccardo/Desktop/codice pacs/build" "/home/riccardo/Desktop/codice pacs/build/tests/example_tests" "/home/riccardo/Desktop/codice pacs/build/tests/example_tests/CMakeFiles/runSimpleTest.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : tests/example_tests/CMakeFiles/runSimpleTest.dir/depend
