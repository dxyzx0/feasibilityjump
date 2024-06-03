cmake_minimum_required(VERSION 3.28)
project(FeasibilityJump_heauristic)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 99)

set(FJ_SOURCES
        src/pbo_fj.cpp
        src/feasibilityjump.cpp
        src/parser/DefaultCallback.cpp
        src/parser/AbcCallback.cpp
)

# Create the C++ library
add_library(
        pbofj SHARED
        ${FJ_SOURCES}
)

# Create the C++ executable
add_executable(
        pbofj_cpp
        main.cpp
)
# Link the C++ library to the C executable
target_link_libraries(pbofj_cpp pbofj)


# Create the C executable
add_executable(
        pbofj_c
        main.c
)
# Link the C++ library to the C executable
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
