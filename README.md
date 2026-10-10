# QPET

QoI-preserving error-bounded lossy compression, as a standalone layer on top of existing compressors.

A lossy compressor bounds the error of the data. An analysis usually looks at a quantity computed from the data, a
quantity of interest (QoI), and a bound on the data says little about the error of a non-linear QoI. QPET
([VLDB'25 paper](https://doi.org/10.14778/3742728.3742739)) compresses so that the QoI itself stays within a tolerance.

This repository is a rewrite of the [paper artifact](https://github.com/JLiu-1/QPET-Artifact). QPET no longer lives
inside a fork of each compressor: it drives unmodified compressors through their public APIs, checks every
decompressed value, and stores the few values that break a bound losslessly.

**Full documentation: [docs/](docs/README.md).**

## Features

- **Backends**: [SZ3](https://github.com/szcompressor/SZ3), [SZo](https://github.com/BingluCS/SZo),
  [ZFP](https://github.com/LLNL/zfp) and [SPERR](https://github.com/NCAR/SPERR), fetched at their latest version and
  built unmodified. Your own compressor can be added with one adapter file.
- **QoI functions** given as text: `x^2`, `log2(x)`, `sin(10x)`, `tanh(x)`, `1/(x+1)`, ... No symbolic algebra
  library is needed.
- **Point-wise QoIs**: the error of `f(x)` at every value is bounded.
- **Regional QoIs**: the error of the mean of `f(x)` over every block of `r` values per dimension is bounded.
- **Optional bound on the data** itself, together with the QoI tolerance.
- **Guarantees that do not depend on the backend**: the bounds are enforced by checking the decompressed data.
- **Automatic choice of the backend's error bound**, with sampling adapted to each compressor and every parameter
  open to the user.
- `float32` and `float64` arrays with 1 to 4 dimensions.
- **Interfaces**: command-line tool, C++, C and Python.

Not supported yet: QoIs of several fields (vector QoIs), and regions other than aligned blocks.

## Build

Needs CMake 3.19 or newer, a C++17 compiler with OpenMP, libzstd (found through pkg-config), and access to GitHub.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
```

The four compressors are downloaded from their official repositories and built into `build/deps`; nothing has to be
prepared by hand. To leave one out, pin a version, or use a local checkout:

```bash
cmake -S . -B build -DQPET_WITH_SZO=OFF -DQPET_ZFP_GIT_TAG=1.0.1 -DQPET_SZ3_SOURCE_DIR=/path/to/SZ3
```

The build produces the tool `build/qpet`, the C++ library `build/libqpet.a`, the C library `build/libqpet_c.so` and
the Python package in `build/python`. See [docs/build.md](docs/build.md).

## Examples

Command line. Dimensions are given fastest-varying first, as in SZ3: a 100x500x500 C array is `-3 500 500 100`.

```bash
# keep x^2 within 0.1 % of its value range at every value; compress with SZ3, then report the errors
build/qpet -b sz3 -f -3 500 500 100 -i Uf48.bin.f32 -q "x^2" -m REL 1e-3 -z U.qpet -a

# keep the mean of x^2 over every 4x4x4 block within tolerance instead, with SPERR
build/qpet -b sperr -f -3 500 500 100 -i Uf48.bin.f32 -q "x^2" -r 4 -m REL 1e-3 -z U.qpet -a

# decompress: the stream records the backend, the data type and the dimensions
build/qpet -z U.qpet -o U.out.f32
```

Python (`export PYTHONPATH=build/python`, or `pip install build/python`):

```python
import numpy as np, qpet

data = np.fromfile("Uf48.bin.f32", dtype=np.float32).reshape(100, 500, 500)
blob = qpet.compress(data, qoi="x^2", qoi_tol=1e-3, backend="sz3")
back = qpet.decompress(blob)
print(len(blob), qpet.verify(data, back, qoi="x^2"))
```

C++ (link the CMake target `qpet`):

```cpp
#include "qpet/qpet.hpp"

qpet::Config conf;
conf.backend = "sperr";
conf.qoi = "log2(x)";
conf.qoi_tol = 1e-3;                                                      // relative to the range of the QoI
std::vector<uint8_t> cmp = qpet::compress(data, {nz, ny, nx}, conf);      // float* or double*, C order

qpet::StreamInfo s = qpet::inspect(cmp.data(), cmp.size());
std::vector<float> dec(s.num);
qpet::decompress(cmp.data(), cmp.size(), dec.data());
```

## Documentation

[Concepts](docs/concepts.md) ·
[Building](docs/build.md) ·
[Command-line tool](docs/cli.md) ·
[QoI expressions](docs/qoi.md) ·
[C++ API](docs/cpp-api.md) ·
[C API](docs/c-api.md) ·
[Python API](docs/python.md) ·
[Tuning](docs/tuning.md) ·
[Backends and adding your own](docs/backends.md) ·
[Stream format](docs/format.md)
