cmake_minimum_required(VERSION 3.28)
project(FeasibilityJump_heauristic)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 99)

set(FJ_SOURCES
        pbo_fj.cpp
        feasibilityjump.cpp
        parser/DefaultCallback.cpp
        parser/AbcCallback.cpp
)

# Create the C++ library
add_library(
        pbo_fj_lib SHARED
        ${FJ_SOURCES}
)

# Create the C++ executable
add_executable(
        pbo_fj
        main.cpp
)
# Link the C++ library to the C executable
target_link_libraries(pbo_fj pbo_fj_lib)


# Create the C executable
add_executable(
        pbo_fj_c
        main.c
)
# Link the C++ library to the C executable
target_link_libraries(pbo_fj_c pbo_fj_lib)


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
    target_compile_definitions(pbo_fj PRIVATE useGMP)
endif ()

target_link_libraries(pbo_fj_lib ${LIBS})
