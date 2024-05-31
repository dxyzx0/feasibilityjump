cmake_minimum_required(VERSION 3.28)
project(FeasibilityJump_heauristic)

set(CMAKE_CXX_STANDARD 23)

# add -DuseGMP to compile with GMP
# add -DuseGMPXX to compile with GMPXX
add_executable(
        pbo_fj
        feasibilityjump.hh
        parser/SimpleParser.h
        parser/DefaultCallback.cpp
        parser/DefaultCallback.h
        parser/AbcCallback.cpp
        parser/AbcCallback.h
        type.h
        pbo_fj.cc
)


#find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS /usr/local/include /usr/local/opt/gmp/include)
find_library(GMP_LIBRARY NAMES gmp PATHS /usr/local/lib /usr/local/opt/gmp/lib)
find_library(GMPXX_LIBRARY NAMES gmpxx PATHS /usr/local/lib /usr/local/opt/gmp/lib)
if (GMP_LIBRARY)
    message(STATUS "GMP Library found: ${GMP_LIBRARY}")
else()
    message(FATAL_ERROR "GMP Library not found")
endif()

if (GMPXX_LIBRARY)
    message(STATUS "GMPXX Library found: ${GMPXX_LIBRARY}")
else()
    message(FATAL_ERROR "GMPXX Library not found")
endif()
include_directories(/usr/local/include)
target_link_libraries(pbo_fj ${GMP_LIBRARY} ${GMPXX_LIBRARY})
target_compile_definitions(pbo_fj PRIVATE useGMP useEigen)

# add Eigen as a dependency
#option(useEigen "Use Eigen" ON)
#if (useEigen)
#    add_executable(pbo_fj
#            feasibilityjump.hh
#            parser/AbcCallback.cpp
#            parser/AbcCallback.h
#            type.h
#
#    )
#    target_link_libraries(pbo_fj ${GMP_LIBRARY} ${GMPXX_LIBRARY})
#    target_compile_definitions(pbo_fj PRIVATE useGMP useEigen)
#endif ()
