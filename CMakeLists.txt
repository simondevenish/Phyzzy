cmake_minimum_required(VERSION 3.15)
project(Phyzzy LANGUAGES C CXX)

# Set C and C++ standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the include directories
include_directories(include)

# Collect source files
set(PHYZZY_SOURCES
    src/physics.c
    src/physicsbody_physics.c
    src/physicsbody.cpp
)

set(PHYZZY_HEADERS
    include/phyzzy/physics.h
    include/phyzzy/physicsbody_physics.h
    include/phyzzy/physicsbody.h
    include/math_types.h
    include/physics_types.h
    include/types.h
    include/vector_math.h
)

# Create a static library
add_library(phyzzy STATIC ${PHYZZY_SOURCES} ${PHYZZY_HEADERS})

# Specify the output directory for the compiled library
set_target_properties(phyzzy PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PhyzzyLib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PhyzzyLib
)

# Add warnings and additional compiler flags for safety
if (CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(phyzzy PRIVATE -Wall -Wextra -Werror -pedantic)
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    target_compile_options(phyzzy PRIVATE /W4 /WX)
endif()

# Optionally, add an executable for testing
option(BUILD_TESTS "Build test executable" OFF)

if (BUILD_TESTS)
    add_executable(test_phyzzy test/test_main.c)
    target_link_libraries(test_phyzzy PRIVATE phyzzy)
    target_include_directories(test_phyzzy PRIVATE include)
endif()
