cmake_minimum_required(VERSION 3.25)
project(FeasibilityJump_heuristic)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 99)

set(FJ_SOURCES
        src/fj_pbo.cpp
        src/fj_solver.cpp
        src/parser/DefaultCallback.cpp
        src/parser/PboCallback.cpp
)

# Create the C++ library
add_library(
        pbofj SHARED
        ${FJ_SOURCES}
)

# Create the C++ executable
add_executable(
        pbofj_cpp
        src/main.cpp
)
target_link_libraries(pbofj_cpp pbofj)


# Create the C executable
add_executable(
        pbofj_c
        src/main.c
)
target_link_libraries(pbofj_c pbofj)


option(USE_GMP "Use GMP library" OFF)

if (USE_GMP)
    include_directories(/usr/local/include)

    find_library(GMP_LIBRARY NAMES gmp PATHS /usr/local/lib /usr/local/opt/gmp/lib)
    find_library(GMPXX_LIBRARY NAMES gmpxx PATHS /usr/local/lib /usr/local/opt/gmp/lib)

    if (GMP_LIBRARY)
        message(STATUS "GMP Library found: ${GMP_LIBRARY}")
    else ()
        message(FATAL_ERROR "GMP Library not found")
    endif ()

    if (GMPXX_LIBRARY)
        message(STATUS "GMPXX Library found: ${GMPXX_LIBRARY}")
    else ()
        message(FATAL_ERROR "GMPXX Library not found")
    endif ()
    set(LIBS ${LIBS} ${GMP_LIBRARY} ${GMPXX_LIBRARY})
    target_compile_definitions(pbofj PRIVATE useGMP)
endif ()

target_link_libraries(pbofj ${LIBS})
