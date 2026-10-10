# Building

## Requirements

- CMake 3.19 or newer
- a C and a C++17 compiler with OpenMP
- libzstd with its headers, found through pkg-config
- git and access to GitHub, unless every compressor is taken from a local checkout
- for the Python package: Python 3.8 or newer and NumPy (not needed to build)

The build was developed and tested on Linux with GCC 13. Other platforms are untested.

## Build and test

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
```

The default build type is `Release`.

## How the compressors are obtained

QPET does not contain SZ3, SZo, ZFP or SPERR. At configure time each one is fetched from its official repository, at
the head of its default branch; at build time it is built unmodified with its own build system and installed into
`build/deps`. Nothing is installed outside the build directory.

| compressor | repository | default branch |
|---|---|---|
| SZ3 | https://github.com/szcompressor/SZ3 | `master` |
| SZo | https://github.com/BingluCS/SZo | `main` |
| ZFP | https://github.com/LLNL/zfp | `develop` |
| SPERR | https://github.com/NCAR/SPERR | `main` |

Because the head of a branch moves, two builds made at different times can contain different versions of a
compressor. `qpet -l` prints the versions of a build, and every compressed stream records the version it was made
with. To pin a version or to build offline:

```bash
# a branch, tag or commit
cmake -S . -B build -DQPET_ZFP_GIT_TAG=1.0.1

# a local checkout (no download)
cmake -S . -B build -DQPET_SZ3_SOURCE_DIR=/path/to/SZ3 -DQPET_SPERR_SOURCE_DIR=/path/to/SPERR
```

## CMake options

| option | default | meaning |
|---|---|---|
| `QPET_WITH_SZ3`, `QPET_WITH_SZO`, `QPET_WITH_ZFP`, `QPET_WITH_SPERR` | `ON` | build the backend |
| `QPET_<NAME>_SOURCE_DIR` | empty | local source tree of the compressor (`<NAME>` is `SZ3`, `SZO`, `ZFP` or `SPERR`) |
| `QPET_<NAME>_GIT_TAG` | default branch | branch, tag or commit to fetch |
| `QPET_SZO_SIMD` | `auto` | SIMD code paths of SZo: `auto` (AVX2 on x86-64), `avx2`, `sve2`, `off`. SZo's streams are the same with and without them |
| `QPET_BUILD_C_API` | `ON` | build `libqpet_c` and the Python package |
| `QPET_BUILD_TOOLS` | `ON` | build the `qpet` command-line tool |
| `QPET_BUILD_TESTS` | `ON` | build the tests |
| `CMAKE_BUILD_TYPE` | `Release` | |

## What the build produces

| path under `build/` | what |
|---|---|
| `qpet` | the [command-line tool](cli.md) |
| `libqpet.a` | the C++ library ([C++ API](cpp-api.md)) |
| `libqpet_c.so` | the [C API](c-api.md), a shared library that contains QPET and the compressors |
| `python/` | the [Python package](python.md) with `libqpet_c.so` next to it |
| `deps/` | the compressors' headers and libraries |
| `test_expr`, `test_roundtrip` | the tests |

There is no `install` target yet.

## Tests

- `expr`: the QoI expression parser, its values and derivatives against reference functions and finite differences,
  and the point-wise error bound estimate.
- `roundtrip`: compress, decompress and verify on synthetic data, for every backend of the build, both data types, 1
  to 4 dimensions, several QoIs, point-wise and regional, including data where the QoI is not finite. The regional
  errors are also checked by code that shares nothing with the library's region handling.
- `python`: the same kind of round trip through the Python package. It is registered only if the Python that CMake
  finds can import NumPy.

## Using QPET from another CMake project

```cmake
add_subdirectory(path/to/QPET qpet)
target_link_libraries(my_cxx_program PRIVATE qpet)     # C++ API
target_link_libraries(my_c_program PRIVATE qpet_c)     # C API
```

Set `QPET_BUILD_TESTS` and `QPET_BUILD_TOOLS` to `OFF` before `add_subdirectory` if you do not want them. The
compressors are built under the top-level build directory, in `deps/`.
