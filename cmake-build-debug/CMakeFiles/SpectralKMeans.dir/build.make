# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.17

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

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2020.3.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2020.3.2\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Tomer\CLionProjects\SpectralKMeans

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Tomer\CLionProjects\SpectralKMeans\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/SpectralKMeans.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SpectralKMeans.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SpectralKMeans.dir/flags.make

CMakeFiles/SpectralKMeans.dir/main.c.obj: CMakeFiles/SpectralKMeans.dir/flags.make
CMakeFiles/SpectralKMeans.dir/main.c.obj: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Tomer\CLionProjects\SpectralKMeans\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/SpectralKMeans.dir/main.c.obj"
	C:\PROGRA~1\MINGW-~1\X86_64~1.0-P\mingw64\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\SpectralKMeans.dir\main.c.obj   -c C:\Users\Tomer\CLionProjects\SpectralKMeans\main.c

CMakeFiles/SpectralKMeans.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SpectralKMeans.dir/main.c.i"
	C:\PROGRA~1\MINGW-~1\X86_64~1.0-P\mingw64\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:\Users\Tomer\CLionProjects\SpectralKMeans\main.c > CMakeFiles\SpectralKMeans.dir\main.c.i

CMakeFiles/SpectralKMeans.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SpectralKMeans.dir/main.c.s"
	C:\PROGRA~1\MINGW-~1\X86_64~1.0-P\mingw64\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:\Users\Tomer\CLionProjects\SpectralKMeans\main.c -o CMakeFiles\SpectralKMeans.dir\main.c.s

# Object files for target SpectralKMeans
SpectralKMeans_OBJECTS = \
"CMakeFiles/SpectralKMeans.dir/main.c.obj"

# External object files for target SpectralKMeans
SpectralKMeans_EXTERNAL_OBJECTS =

SpectralKMeans.exe: CMakeFiles/SpectralKMeans.dir/main.c.obj
SpectralKMeans.exe: CMakeFiles/SpectralKMeans.dir/build.make
SpectralKMeans.exe: CMakeFiles/SpectralKMeans.dir/linklibs.rsp
SpectralKMeans.exe: CMakeFiles/SpectralKMeans.dir/objects1.rsp
SpectralKMeans.exe: CMakeFiles/SpectralKMeans.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Tomer\CLionProjects\SpectralKMeans\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable SpectralKMeans.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\SpectralKMeans.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SpectralKMeans.dir/build: SpectralKMeans.exe

.PHONY : CMakeFiles/SpectralKMeans.dir/build

CMakeFiles/SpectralKMeans.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\SpectralKMeans.dir\cmake_clean.cmake
.PHONY : CMakeFiles/SpectralKMeans.dir/clean

CMakeFiles/SpectralKMeans.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Tomer\CLionProjects\SpectralKMeans C:\Users\Tomer\CLionProjects\SpectralKMeans C:\Users\Tomer\CLionProjects\SpectralKMeans\cmake-build-debug C:\Users\Tomer\CLionProjects\SpectralKMeans\cmake-build-debug C:\Users\Tomer\CLionProjects\SpectralKMeans\cmake-build-debug\CMakeFiles\SpectralKMeans.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SpectralKMeans.dir/depend

