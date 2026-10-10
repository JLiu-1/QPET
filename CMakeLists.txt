cmake_minimum_required(VERSION 3.19)
project(QPET VERSION 2.0.0 LANGUAGES C CXX)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type to build" FORCE)
endif ()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(QPET_WITH_SZ3 "build the SZ3 backend" ON)
option(QPET_WITH_SZO "build the SZo backend" ON)
option(QPET_WITH_ZFP "build the ZFP backend" ON)
option(QPET_WITH_SPERR "build the SPERR backend" ON)
option(QPET_BUILD_C_API "build the C interface (libqpet_c) and the Python package" ON)
option(QPET_BUILD_TOOLS "build the qpet command-line tool" ON)
option(QPET_BUILD_TESTS "build the tests" ON)

find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_search_module(ZSTD REQUIRED IMPORTED_TARGET libzstd)

include(cmake/QPETDeps.cmake)

add_library(qpet STATIC src/expr.cpp src/core.cpp src/registry.cpp)
add_library(QPET::qpet ALIAS qpet)
target_include_directories(qpet PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(qpet PRIVATE PkgConfig::ZSTD OpenMP::OpenMP_CXX)
target_compile_options(qpet PRIVATE -Wall -Wextra)

# One translation unit per backend: the compressors' headers do not all coexist (SZ3 and SZo both define
# versionInt() and versionStr() in the global namespace).
foreach (backend SZ3 SZO ZFP SPERR)
    if (QPET_WITH_${backend})
        string(TOLOWER ${backend} name)
        target_sources(qpet PRIVATE src/backends/${name}.cpp)
        target_compile_definitions(qpet PRIVATE QPET_WITH_${backend})
        target_link_libraries(qpet PRIVATE qpet_dep_${name})
        if (TARGET qpet_ext_${name})  # the compressor is built by this build, see cmake/QPETDeps.cmake
            add_dependencies(qpet qpet_ext_${name})
        endif ()
    endif ()
endforeach ()

if (QPET_WITH_SZO)
    # SZo is not header-only: its entropy coder needs the FSE sources it vendors. Hidden visibility keeps their FSE_*
    # symbols apart from the ones inside libzstd.
    target_sources(qpet PRIVATE ${QPET_SZO_ZFSE_SOURCES})
    set_source_files_properties(${QPET_SZO_ZFSE_SOURCES} PROPERTIES COMPILE_OPTIONS "-fvisibility=hidden;-w")
    if (QPET_SZO_SIMD_FLAGS)
        set_source_files_properties(src/backends/szo.cpp PROPERTIES COMPILE_OPTIONS "${QPET_SZO_SIMD_FLAGS}")
    endif ()
endif ()

if (QPET_BUILD_C_API)
    # The C interface as a shared library that holds QPET and the compressors. The Python package loads it.
    add_library(qpet_c SHARED src/capi.cpp)
    target_link_libraries(qpet_c PRIVATE qpet)
    target_include_directories(qpet_c PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

    # <build>/python is the Python package with the library next to it: use it through PYTHONPATH or pip install it.
    set(QPET_PYTHON_DIR ${CMAKE_BINARY_DIR}/python)
    add_custom_command(TARGET qpet_c POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory ${QPET_PYTHON_DIR}/qpet
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/pyproject.toml ${QPET_PYTHON_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/qpet/__init__.py ${QPET_PYTHON_DIR}/qpet
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:qpet_c> ${QPET_PYTHON_DIR}/qpet)
endif ()

if (QPET_BUILD_TOOLS)
    add_executable(qpet_cli tools/qpet.cpp)
    set_target_properties(qpet_cli PROPERTIES OUTPUT_NAME qpet)
    target_link_libraries(qpet_cli PRIVATE qpet)
endif ()

if (QPET_BUILD_TESTS)
    enable_testing()
    add_executable(test_expr tests/test_expr.cpp)
    target_link_libraries(test_expr PRIVATE qpet)
    add_test(NAME expr COMMAND test_expr)
    add_executable(test_roundtrip tests/test_roundtrip.cpp)
    target_link_libraries(test_roundtrip PRIVATE qpet)
    add_test(NAME roundtrip COMMAND test_roundtrip)
    if (QPET_BUILD_C_API)
        find_package(Python3 COMPONENTS Interpreter)
        if (Python3_FOUND)
            execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy" RESULT_VARIABLE QPET_NO_NUMPY
                    OUTPUT_QUIET ERROR_QUIET)
            if (NOT QPET_NO_NUMPY)
                add_test(NAME python COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_python.py)
                set_tests_properties(python PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
            else ()
                message(STATUS "QPET: numpy not found, the Python test is not registered")
            endif ()
        endif ()
    endif ()
endif ()
