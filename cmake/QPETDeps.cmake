# The compressors QPET wraps. Each is fetched from its official repository at the head of its default branch (or taken
# from a local checkout), built unmodified with its own build system, and installed into a private prefix.
#
#   -DQPET_<NAME>_SOURCE_DIR=/path   use a local checkout instead of fetching
#   -DQPET_<NAME>_GIT_TAG=<ref>      fetch this branch, tag or commit

include(FetchContent)
include(ExternalProject)

set(QPET_DEPS_PREFIX ${CMAKE_BINARY_DIR}/deps)
file(MAKE_DIRECTORY ${QPET_DEPS_PREFIX}/include)

set(QPET_SZ3_GIT_REPOSITORY https://github.com/szcompressor/SZ3.git)
set(QPET_SZO_GIT_REPOSITORY https://github.com/BingluCS/SZo.git)
set(QPET_ZFP_GIT_REPOSITORY https://github.com/LLNL/zfp.git)
set(QPET_SPERR_GIT_REPOSITORY https://github.com/NCAR/SPERR.git)
set(QPET_SZ3_DEFAULT_TAG master)
set(QPET_SZO_DEFAULT_TAG main)
set(QPET_ZFP_DEFAULT_TAG develop)
set(QPET_SPERR_DEFAULT_TAG main)

# Sets QPET_<NAME>_SRC to the compressor's source tree, downloading it if needed. The tree is not added to this build.
function(qpet_get_source NAME)
    set(QPET_${NAME}_SOURCE_DIR "" CACHE PATH "local ${NAME} checkout to build instead of fetching")
    set(QPET_${NAME}_GIT_TAG ${QPET_${NAME}_DEFAULT_TAG} CACHE STRING "${NAME} branch, tag or commit to fetch")
    if (QPET_${NAME}_SOURCE_DIR)
        get_filename_component(src "${QPET_${NAME}_SOURCE_DIR}" ABSOLUTE)
        if (NOT EXISTS ${src}/CMakeLists.txt)
            message(FATAL_ERROR "QPET_${NAME}_SOURCE_DIR=${src} is not a source tree")
        endif ()
    else ()
        string(TOLOWER ${NAME} lower)
        FetchContent_Declare(qpet_src_${lower}
                GIT_REPOSITORY ${QPET_${NAME}_GIT_REPOSITORY}
                GIT_TAG ${QPET_${NAME}_GIT_TAG}
                SOURCE_SUBDIR qpet-does-not-add-this-to-the-build)
        FetchContent_MakeAvailable(qpet_src_${lower})
        set(src ${qpet_src_${lower}_SOURCE_DIR})
    endif ()
    message(STATUS "QPET: ${NAME} source: ${src}")
    set(QPET_${NAME}_SRC ${src} PARENT_SCOPE)
endfunction()

function(qpet_build NAME)
    string(TOLOWER ${NAME} lower)
    ExternalProject_Add(qpet_ext_${lower}
            SOURCE_DIR ${QPET_${NAME}_SRC}
            PREFIX ${CMAKE_BINARY_DIR}/ext/${lower}
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${QPET_DEPS_PREFIX}
            -DCMAKE_INSTALL_LIBDIR=lib
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DBUILD_SHARED_LIBS=OFF
            -DBUILD_TESTING=OFF
            ${ARGN})
endfunction()

if (QPET_WITH_SZ3)
    qpet_get_source(SZ3)
    qpet_build(SZ3 -DBUILD_SZ3_BINARY=OFF)
    add_library(qpet_dep_sz3 INTERFACE)
    target_include_directories(qpet_dep_sz3 INTERFACE ${QPET_DEPS_PREFIX}/include)
    target_link_libraries(qpet_dep_sz3 INTERFACE PkgConfig::ZSTD OpenMP::OpenMP_CXX)
endif ()

if (QPET_WITH_SZO)
    qpet_get_source(SZO)
    qpet_build(SZO -DBUILD_SZo_BINARY=OFF -DBUILD_SZOC_LIBRARY=OFF -DSZo_USE_BUNDLED_ZSTD=OFF)
    add_library(qpet_dep_szo INTERFACE)
    target_include_directories(qpet_dep_szo INTERFACE ${QPET_DEPS_PREFIX}/include)
    target_link_libraries(qpet_dep_szo INTERFACE PkgConfig::ZSTD OpenMP::OpenMP_CXX)
    file(GLOB QPET_SZO_ZFSE_SOURCES ${QPET_SZO_SRC}/include/SZo/encoder/zfse/*.c)
    # SZo's SIMD paths are chosen by the flags of the translation unit that includes it. Its streams are the same
    # with and without them.
    set(QPET_SZO_SIMD "auto" CACHE STRING "SIMD for the SZo backend: auto, avx2, sve2 or off")
    set(QPET_SZO_SIMD_FLAGS "")
    if (QPET_SZO_SIMD STREQUAL "avx2" OR (QPET_SZO_SIMD STREQUAL "auto" AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64"))
        set(QPET_SZO_SIMD_FLAGS "-mavx2;-mfma")
    elseif (QPET_SZO_SIMD STREQUAL "sve2")
        set(QPET_SZO_SIMD_FLAGS "-march=armv8.6-a+sve2")
    endif ()
endif ()

if (QPET_WITH_ZFP)
    qpet_get_source(ZFP)
    qpet_build(ZFP -DBUILD_UTILITIES=OFF -DZFP_WITH_OPENMP=ON
            BUILD_BYPRODUCTS ${QPET_DEPS_PREFIX}/lib/libzfp.a)
    add_library(qpet_dep_zfp STATIC IMPORTED)
    set_target_properties(qpet_dep_zfp PROPERTIES
            IMPORTED_LOCATION ${QPET_DEPS_PREFIX}/lib/libzfp.a
            INTERFACE_INCLUDE_DIRECTORIES ${QPET_DEPS_PREFIX}/include
            INTERFACE_LINK_LIBRARIES "OpenMP::OpenMP_C;m")
endif ()

if (QPET_WITH_SPERR)
    qpet_get_source(SPERR)
    qpet_build(SPERR -DBUILD_UNIT_TESTS=OFF -DBUILD_CLI_UTILITIES=OFF -DUSE_OMP=ON
            BUILD_BYPRODUCTS ${QPET_DEPS_PREFIX}/lib/libSPERR.a)
    add_library(qpet_dep_sperr STATIC IMPORTED)
    set_target_properties(qpet_dep_sperr PROPERTIES
            IMPORTED_LOCATION ${QPET_DEPS_PREFIX}/lib/libSPERR.a
            INTERFACE_INCLUDE_DIRECTORIES ${QPET_DEPS_PREFIX}/include
            INTERFACE_LINK_LIBRARIES OpenMP::OpenMP_CXX)
endif ()
