# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.20

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2021.2.3\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2021.2.3\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\gemma\CLionProjects\metodesnumericsII

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\gemma\CLionProjects\metodesnumericsII\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/metodesnumericsII.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/metodesnumericsII.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/metodesnumericsII.dir/flags.make

CMakeFiles/metodesnumericsII.dir/main.c.obj: CMakeFiles/metodesnumericsII.dir/flags.make
CMakeFiles/metodesnumericsII.dir/main.c.obj: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\gemma\CLionProjects\metodesnumericsII\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/metodesnumericsII.dir/main.c.obj"
	C:\MinGW\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles\metodesnumericsII.dir\main.c.obj -c C:\Users\gemma\CLionProjects\metodesnumericsII\main.c

CMakeFiles/metodesnumericsII.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/metodesnumericsII.dir/main.c.i"
	C:\MinGW\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:\Users\gemma\CLionProjects\metodesnumericsII\main.c > CMakeFiles\metodesnumericsII.dir\main.c.i

CMakeFiles/metodesnumericsII.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/metodesnumericsII.dir/main.c.s"
	C:\MinGW\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:\Users\gemma\CLionProjects\metodesnumericsII\main.c -o CMakeFiles\metodesnumericsII.dir\main.c.s

# Object files for target metodesnumericsII
metodesnumericsII_OBJECTS = \
"CMakeFiles/metodesnumericsII.dir/main.c.obj"

# External object files for target metodesnumericsII
metodesnumericsII_EXTERNAL_OBJECTS =

metodesnumericsII.exe: CMakeFiles/metodesnumericsII.dir/main.c.obj
metodesnumericsII.exe: CMakeFiles/metodesnumericsII.dir/build.make
metodesnumericsII.exe: CMakeFiles/metodesnumericsII.dir/linklibs.rsp
metodesnumericsII.exe: CMakeFiles/metodesnumericsII.dir/objects1.rsp
metodesnumericsII.exe: CMakeFiles/metodesnumericsII.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\gemma\CLionProjects\metodesnumericsII\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable metodesnumericsII.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\metodesnumericsII.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/metodesnumericsII.dir/build: metodesnumericsII.exe
.PHONY : CMakeFiles/metodesnumericsII.dir/build

CMakeFiles/metodesnumericsII.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\metodesnumericsII.dir\cmake_clean.cmake
.PHONY : CMakeFiles/metodesnumericsII.dir/clean

CMakeFiles/metodesnumericsII.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\gemma\CLionProjects\metodesnumericsII C:\Users\gemma\CLionProjects\metodesnumericsII C:\Users\gemma\CLionProjects\metodesnumericsII\cmake-build-debug C:\Users\gemma\CLionProjects\metodesnumericsII\cmake-build-debug C:\Users\gemma\CLionProjects\metodesnumericsII\cmake-build-debug\CMakeFiles\metodesnumericsII.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/metodesnumericsII.dir/depend

