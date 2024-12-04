cmake_minimum_required(VERSION 3.15)
project(Phyzzy LANGUAGES C CXX)

# Set C and C++ standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the include directories
# Add both the root "include" directory and "include/phyzzy" for simplicity
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/phyzzy)

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
add_library(PhyzzyLib STATIC ${PHYZZY_SOURCES} ${PHYZZY_HEADERS})

# Specify the output directory for the compiled library
set_target_properties(PhyzzyLib PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PhyzzyLib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PhyzzyLib
)

# Add warnings and additional compiler flags for safety
if (CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(PhyzzyLib PRIVATE -Wall -Wextra -Werror -pedantic)
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    target_compile_options(PhyzzyLib PRIVATE /W4 /WX)
endif()

# Optionally, add an executable for testing
option(BUILD_TESTS "Build test executable" OFF)

if (BUILD_TESTS)
    add_executable(test_phyzzy test/test_main.c)
    target_link_libraries(test_phyzzy PRIVATE PhyzzyLib)
    target_include_directories(test_phyzzy PRIVATE ${CMAKE_SOURCE_DIR}/include)
endif()
