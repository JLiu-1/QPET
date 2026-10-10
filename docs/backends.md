# Backends

A backend is an existing error-bounded lossy compressor that QPET drives through its public API. QPET gives it one
absolute error bound for the whole array, and does not modify it. `qpet -l` lists the backends of a build with the
versions of their compressors.

| backend | compressor | mode QPET uses | data |
|---|---|---|---|
| `sz3` | [SZ3](https://github.com/szcompressor/SZ3) | absolute error bound | `float32`, `float64`; 1-D to 4-D |
| `szo` | [SZo](https://github.com/BingluCS/SZo) | absolute error bound | `float32`, `float64`; 1-D to 4-D |
| `zfp` | [ZFP](https://github.com/LLNL/zfp) | fixed-accuracy | `float32`, `float64`; 1-D to 4-D |
| `sperr` | [SPERR](https://github.com/NCAR/SPERR) | point-wise error (PWE) | `float32`, `float64`; 2-D and 3-D |

## Notes on each

### `sz3`

- Option `algo`: `interp_lorenzo` (SZ3's default: interpolation or Lorenzo, chosen by SZ3), `interp`, `lorenzo_reg`.
- Single-threaded unless `backend_parallel` (`-P`) is set. SZ3's OpenMP mode splits the array and compresses the
  parts separately, which is faster and lowers the compression ratio.

### `szo`

- Same `algo` option and `backend_parallel` behaviour as `sz3`.
- SZo compresses in place and overwrites its input, so the adapter gives it a copy. This costs one more copy of the
  input in memory during compression.
- Its SIMD code paths are selected at build time (`QPET_SZO_SIMD`, see [Building](build.md)).

### `zfp`

- The error bound is the tolerance of ZFP's fixed-accuracy mode. ZFP only uses the binary exponent of the tolerance,
  and its actual errors are usually far below it, so the tolerance QPET ends up with is often much larger than the
  error bound another backend would get.
- Compression uses ZFP's OpenMP execution with the given threads; the stream does not depend on the thread count
  (checked on one 3-D field). Decompression is serial.
- The tolerance is stored in QPET's header, because ZFP needs it to decode.

### `sperr`

- 2-D and 3-D only.
- Option `chunk`: edge of the chunks a 3-D volume is split into (default 256). Chunks are compressed in parallel
  with the given threads; the stream does not depend on the thread count (checked on one 3-D field).
- 2-D slices are compressed without SPERR's own header; QPET's header carries the dimensions.
- SPERR cannot compress data that holds NaN or infinite values, or magnitudes above about 2e31 (such as a fill value of 1e36); QPET reports an error for such data. The other backends accept it.

## Passing options

| interface | form |
|---|---|
| command line | `-O algo=interp` |
| C++ | `conf.backend_options["algo"] = "interp";` |
| C | `backend.algo=interp` in the options text |
| Python | `backend_options={"algo": "interp"}` |

An option that a backend does not know is ignored.

## Adding your own compressor

A compressor is connected to QPET by a backend adapter: one source file that calls the compressor's API. Everything
else (tuning, validation, the stream format, the command-line tool, the benchmarks) then works with it. The compressor
itself is not modified. It needs to

- take one absolute error bound, for `float` and/or `double` arrays;
- be callable from several threads at once (the tuning compresses sample blocks in parallel);
- leave its input unchanged (if it does not, copy the input in the adapter, as `src/backends/szo.cpp` does).

It does not need to respect the error bound strictly: QPET checks every decompressed value.

The steps below add a compressor called `mycomp`. `src/backends/zfp.cpp` is a complete adapter of about 100 lines to
start from.

**1. Write `src/backends/mycomp.cpp`**

```cpp
#include "mycomp.h"            // the compressor's own header
#include "qpet/backend.hpp"

namespace qpet {

namespace {

class MyCompBackend : public Backend {
   public:
    std::string name() const override { return "mycomp"; }      // the name given to -b
    std::string version() const override { return MYCOMP_VERSION; }   // stored in every stream

    bool supports(DType dtype, size_t ndims) const override { return dtype == DType::F32 && ndims <= 3; }

    // dims are in C order (slowest-varying first). params.abs_eb is the error bound, params.threads the threads the
    // compressor may use, params.options the -O key=value pairs of the user.
    std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims,
                                  const BackendParams &params) const override {
        std::vector<uint8_t> out(mycomp_bound(dims));
        out.resize(mycomp_compress(static_cast<const float *>(data), dims, params.abs_eb, out.data()));
        return out;
    }

    // out holds prod(dims) values. dtype, dims and params.abs_eb are the ones compress() was called with: QPET
    // stores them, so the compressor's stream does not have to.
    void decompress(const uint8_t *buf, size_t size, DType dtype, const std::vector<size_t> &dims,
                    const BackendParams &params, void *out) const override {
        mycomp_decompress(buf, size, dims, static_cast<float *>(out));
    }
};

}  // namespace

const Backend &mycomp_backend() {
    static const MyCompBackend b;
    return b;
}

}  // namespace qpet
```

Report a failure by throwing a `std::exception`.

**2. Register it in `src/registry.cpp`**

```cpp
const Backend &mycomp_backend();          // next to the other declarations
...
#ifdef QPET_WITH_MYCOMP                   // in backends()
        &mycomp_backend(),
#endif
```

**3. Add it to the build, in `CMakeLists.txt`**

```cmake
option(QPET_WITH_MYCOMP "build the mycomp backend" ON)
...
foreach (backend SZ3 SZO ZFP SPERR MYCOMP)
```

The loop compiles `src/backends/mycomp.cpp`, defines `QPET_WITH_MYCOMP` and links the target `qpet_dep_mycomp`.

**4. Define `qpet_dep_mycomp`, in `cmake/QPETDeps.cmake`**

It is the target that carries the compressor's include directories and libraries. For a compressor that is already
installed:

```cmake
if (QPET_WITH_MYCOMP)
    find_package(mycomp REQUIRED)
    add_library(qpet_dep_mycomp INTERFACE)
    target_link_libraries(qpet_dep_mycomp INTERFACE mycomp::mycomp)
endif ()
```

To have QPET fetch and build it instead, as it does for the four built-in ones, follow the `ZFP` block in that file:
`qpet_get_source` downloads the source (or takes `-DQPET_MYCOMP_SOURCE_DIR=...`), `qpet_build` builds and installs it
into `build/deps`, and an imported library points at the result.

**5. Build and try it**

```bash
cmake -S . -B build && cmake --build build -j
build/qpet -l                                    # lists mycomp and its version
build/qpet -b mycomp -f -3 512 512 512 -i data.f32 -q "x^2" -m REL 1e-3 -a -v
ctest --test-dir build                           # the round-trip test runs on every backend of the build
```

`-a` reports the compression ratio, the speeds and whether the bounds hold; `-v` shows what the tuning tried.

**Optional: tuning defaults for your compressor**

The tuning compresses a sample instead of the whole input, and the general defaults assume nothing about the
compressor. If you know how yours behaves, override `tuning_defaults()` in the adapter to return the
[tuning parameters](tuning.md) that suit it; a user can still override each of them. For example, a compressor
that codes 4x4x4 cells independently and reacts only to the binary exponent of the error bound:

```cpp
    TuningParams tuning_defaults(size_t /*ndims*/) const override {
        TuningParams p;
        p.independent_cell = 4;
        p.pow2_error_bounds = 1;
        return p;
    }
```

To judge a choice of parameters, compare what the tuning picks with a scan of the error bound on a few of your
inputs: `qpet -E <bound>` compresses with a fixed one, `qpet -v` shows what the tuning tried.

**Why one file per compressor**

The compressors' headers do not all coexist in one translation unit: SZ3 and SZo both define `versionInt()` and
`versionStr()` in the global namespace.
