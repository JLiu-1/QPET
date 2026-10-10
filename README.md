QPET: Quantity-of-Interest-Preserving Error-Bounded Lossy Compression
=====

This repository implements the method of the paper

> Jinyang Liu, Pu Jiao, Kai Zhao, Xin Liang, Sheng Di, and Franck Cappello.
> **QPET: A Versatile and Portable Quantity-of-Interest-Preservation Framework for Error-Bounded Lossy Compression.**
> PVLDB 18(8): 2440-2453, 2025. [doi:10.14778/3742728.3742739](https://doi.org/10.14778/3742728.3742739)

An error-bounded lossy compressor bounds the error of the data. An analysis usually looks at a quantity computed from
the data, a quantity of interest (QoI), and a bound on the data says little about the error of a non-linear QoI. QPET
compresses floating-point arrays so that the QoI itself stays within the tolerance you set. It works on top of an
existing compressor, [SZ3](https://github.com/szcompressor/SZ3), [SZo](https://github.com/BingluCS/SZo),
[ZFP](https://github.com/LLNL/zfp) or [SPERR](https://github.com/NCAR/SPERR), without modifying it. It can be used
from the `qpet` command line, C++, C and Python.

The code is a rewrite of the [paper's artifact](https://github.com/JLiu-1/QPET-Artifact), made with
[Claude Code](https://claude.com/claude-code) (Claude Opus 5.5) under the direction of the paper's first author. It
differs from the artifact in some respects, listed [below](#differences-from-the-paper).

**Full documentation: [docs/](docs/README.md).**

## Installation

You need CMake 3.19 or newer, a C++17 compiler with OpenMP, libzstd (found through pkg-config), and access to GitHub.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
```

The four compressors are downloaded from their official repositories at their latest version and built into
`build/deps`; nothing has to be prepared by hand. The tool is `build/qpet`, the C++ library `build/libqpet.a`, the C
library `build/libqpet_c.so`, and the Python package is in `build/python`.

To leave a compressor out, pin its version, or use a local checkout of it:

```bash
cmake -S . -B build -DQPET_WITH_SZO=OFF -DQPET_ZFP_GIT_TAG=1.0.1 -DQPET_SZ3_SOURCE_DIR=/path/to/SZ3
```

There is no install step yet. All build options are in [docs/build.md](docs/build.md).

## Quick start

The examples use the field `Uf48.bin.f32` (100x500x500, `float32`) of the Hurricane Isabel dataset from
[SDRBench](https://sdrbench.github.io/).

### Command line

```bash
# Keep x^2 within 0.1 % of its value range at every value, with SZ3. -a decompresses and prints the errors.
build/qpet -b sz3 -f -3 500 500 100 -i Uf48.bin.f32 -q "x^2" -m REL 1e-3 -z U.qpet -a
# Keep the mean of x^2 over every 4x4x4 block within the tolerance instead, with SPERR
build/qpet -b sperr -f -3 500 500 100 -i Uf48.bin.f32 -q "x^2" -r 4 -m REL 1e-3 -z U.qpet -a
# Decompress. The compressed file records the compressor, the data type and the dimensions.
build/qpet -z U.qpet -o U.out.f32
```

Run `qpet` without arguments for every option; the common ones:

| Option | Meaning |
|---|---|
| `-b <backend>` | compressor: `sz3`, `szo`, `zfp`, `sperr` (`qpet -l` lists them with their versions) |
| `-f`, `-d` | data type: float, double |
| `-1 nx`, `-2 nx ny`, `-3 nx ny nz`, `-4 nx ny nz nw` | dimensions, fastest-varying first: `-3 nx ny nz` is `data[nz][ny][nx]` |
| `-i <file>` | original data, raw binary |
| `-z <file>` | compressed file: written when compressing, read when decompressing |
| `-o <file>` | decompressed data, raw binary |
| `-q <QoI>` | the QoI as a function of one value: `"x^2"`, `"log2(x)"`, `"tanh(x)"`, `"1/(x+1)"`, piecewise as `"x > 0 ? log2(x) : 0"`, ... |
| `-r <edge>` | regional QoI: bound the error of the mean of the QoI over blocks of `<edge>` values per dimension |
| `-m <mode> <tol>` | QoI error tolerance: `ABS`, or `REL` (relative to the value range of the QoI) |
| `-M <mode> <bound>` | optional bound on the error of the data itself: `ABS` or `REL` |
| `-a` | after compression, decompress and print error statistics |
| `-t <n>` | threads |

### Python

```python
import numpy as np, qpet        # export PYTHONPATH=build/python, or pip install build/python

data = np.fromfile("Uf48.bin.f32", dtype=np.float32).reshape(100, 500, 500)
blob = qpet.compress(data, qoi="x^2", qoi_tol=1e-3, backend="sz3")    # bytes
back = qpet.decompress(blob)                                          # same shape and dtype as data
print(data.nbytes / len(blob), qpet.verify(data, back, qoi="x^2"))
```

### C++

```cpp
#include <qpet/qpet.hpp>
#include <vector>

int main() {
    std::vector<float> data(100 * 500 * 500, 1.0f);
    qpet::Config conf;
    conf.backend = "sperr";
    conf.qoi = "log2(x)";
    conf.qoi_tol = 1e-3;  // relative to the value range of the QoI, unless qoi_tol_mode is ABS

    // 500 is the fastest-varying dimension
    std::vector<uint8_t> cmp = qpet::compress(data.data(), {100, 500, 500}, conf);

    qpet::StreamInfo info = qpet::inspect(cmp.data(), cmp.size());  // read from the compressed data
    std::vector<float> dec(info.num);
    qpet::decompress(cmp.data(), cmp.size(), dec.data());
}
```

In a CMake project, `add_subdirectory` this repository and link the target `qpet`.

## What QPET guarantees

For the QoI function `f`, every original value `x` and its decompressed value `x'`:

* point-wise QoI (the default): `|f(x) - f(x')| <= tolerance` for every value;
* regional QoI (`-r`): `|mean f(x) - mean f(x')| <= tolerance` over every block;
* with a data bound (`-M`): `|x - x'| <= bound` for every value;
* where `f(x)` is not finite, for example `log2` of a value that is not positive, the value is reproduced exactly.

These hold whichever compressor is used, and even if it does not respect its own error bound: QPET decompresses what
the compressor produced, checks it against the original, and stores the values that break a bound losslessly.
[docs/concepts.md](docs/concepts.md) has the exact definitions.

## Interfaces

| Interface | How to use / where |
|---|---|
| Command line | `build/qpet`; see [above](#command-line) and [docs/cli.md](docs/cli.md) |
| C++ | `#include <qpet/qpet.hpp>`, CMake target `qpet`; [docs/cpp-api.md](docs/cpp-api.md) |
| C | `#include <qpet/qpet_c.h>`, library `qpet_c`; [docs/c-api.md](docs/c-api.md) |
| Python | package `qpet` in `build/python`, works on NumPy arrays; [docs/python.md](docs/python.md) |

## Supported compressors

| Backend | Compressor | Data |
|---|---|---|
| `sz3` | [SZ3](https://github.com/szcompressor/SZ3) | float and double, 1-D to 4-D |
| `szo` | [SZo](https://github.com/BingluCS/SZo) | float and double, 1-D to 4-D |
| `zfp` | [ZFP](https://github.com/LLNL/zfp), fixed-accuracy mode | float and double, 1-D to 4-D |
| `sperr` | [SPERR](https://github.com/NCAR/SPERR), point-wise error mode | float and double, 2-D and 3-D |

Another compressor is added with one adapter file; [docs/backends.md](docs/backends.md) walks through it.

## Differences from the paper

* **One error bound per array.** In the paper, SZ3-QPET and HPEZ-QPET also quantize a separate error bound for each
  value inside the compressor. Here the compressor is unmodified and receives a single error bound, as SPERR-QPET
  does in the paper. Where the tolerance of the values varies strongly in space, the paper's SZ3-QPET compresses
  better.
* **QoIs of several fields are not supported yet** (the paper's vector QoIs, such as a velocity magnitude). Point-wise
  QoIs of one value and their means over aligned blocks are.
* **HPEZ is not a backend.**
* **No SymEngine.** QoI functions are parsed by QPET itself and differentiated automatically.
* **The error bound given to the compressor is chosen differently**: by trying candidates on a sample of the input,
  with sampling adapted to each compressor. [docs/tuning.md](docs/tuning.md) describes it.

The numbers of the paper come from its [artifact](https://github.com/JLiu-1/QPET-Artifact), not from this code.

## Further documentation

The documentation starts at [docs/README.md](docs/README.md), which lists every page:

* [docs/concepts.md](docs/concepts.md): what is guaranteed, point-wise and regional QoIs, tolerances, how the
  compression works, limits.
* [docs/qoi.md](docs/qoi.md): the syntax of QoI functions.
* [docs/cli.md](docs/cli.md), [docs/cpp-api.md](docs/cpp-api.md), [docs/c-api.md](docs/c-api.md),
  [docs/python.md](docs/python.md): the four interfaces in full.
* [docs/tuning.md](docs/tuning.md): how the compressor's error bound is chosen, and the parameters that control it.
* [docs/backends.md](docs/backends.md): notes on each compressor, and how to add your own.
* [docs/build.md](docs/build.md): build options, pinning compressor versions, tests.
* [docs/format.md](docs/format.md): the layout of a compressed stream.

## Citing QPET

```bibtex
@article{liu2025qpet,
  author  = {Liu, Jinyang and Jiao, Pu and Zhao, Kai and Liang, Xin and Di, Sheng and Cappello, Franck},
  title   = {{QPET}: A Versatile and Portable Quantity-of-Interest-Preservation Framework for Error-Bounded Lossy Compression},
  journal = {Proceedings of the VLDB Endowment},
  volume  = {18},
  number  = {8},
  pages   = {2440--2453},
  year    = {2025},
  doi     = {10.14778/3742728.3742739}
}
```

Please also cite the compressor you use QPET with.

## License and contact

QPET is released under the BSD 3-Clause license; see [LICENSE](LICENSE). The compressors it downloads and links keep
their own licenses.

Report bugs and ask questions in [GitHub issues](https://github.com/JLiu-1/QPET/issues).
