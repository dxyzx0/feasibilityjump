cmake_minimum_required(VERSION 3.28)
project(FeasibilityJump_heauristic)

set(CMAKE_CXX_STANDARD 23)

add_executable(
        pbo_fj
        pbo_fj.cpp
        feasibilityjump.h
        type.h
        parser/SimpleParser.h
        parser/DefaultCallback.cpp
        parser/DefaultCallback.h
        parser/AbcCallback.cpp
        parser/AbcCallback.h
)

option(USE_GMP "Use GMP library" OFF)

if (USE_GMP)
    include_directories(/usr/local/include)

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
    set(LIBS ${LIBS} ${GMP_LIBRARY} ${GMPXX_LIBRARY})
    target_compile_definitions(pbo_fj PRIVATE useGMP)
endif ()

target_link_libraries(pbo_fj ${LIBS})
