cmake_minimum_required(VERSION 3.28)
project(FeasibilityJump_heauristic)

set(CMAKE_CXX_STANDARD 23)

# Set the COPT_HOME environment variable
set(COPT_HOME $ENV{COPT_HOME})
message(STATUS "COPT_HOME: ${COPT_HOME}")

# Add the include directories
include_directories(${COPT_HOME}/include)
include_directories(${COPT_HOME}/include/coptcpp_inc)

# Find the library
find_library(COPT_LIBRARY NAMES copt_cpp PATHS ${COPT_HOME}/lib)
if (NOT COPT_LIBRARY)
    message(FATAL_ERROR "COPT library not found")
endif()

add_executable(
        fjcopt
        copt_fj.cc
)

# Link the library
target_link_libraries(fjcopt ${COPT_LIBRARY})
