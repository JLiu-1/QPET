# C++ API

```cpp
#include "qpet/qpet.hpp"
```

Everything is in the namespace `qpet`. Link the CMake target `qpet` (see [Building](build.md)). Functions report
failures by throwing: `std::invalid_argument` for bad arguments or options, `std::runtime_error` for a corrupted
stream or a failure inside a backend.

Arrays are given as a pointer plus dimensions in **C order** (slowest-varying first): a 100x200x300 C array has
`dims = {100, 200, 300}`. There are 1 to 4 dimensions.

## Compress

```cpp
std::vector<uint8_t> compress(const float  *data, const std::vector<size_t> &dims, const Config &conf, CompressInfo *info = nullptr);
std::vector<uint8_t> compress(const double *data, const std::vector<size_t> &dims, const Config &conf, CompressInfo *info = nullptr);
std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims, const Config &conf, CompressInfo *info = nullptr);
```

`DType` is `DType::F32` or `DType::F64`. The input is not modified. The result is the complete stream.

### `Config`

| field | default | meaning |
|---|---|---|
| `std::string backend` | `"sz3"` | `"sz3"`, `"szo"`, `"zfp"`, `"sperr"` |
| `std::string qoi` | `"x"` | the QoI function, see [QoI expressions](qoi.md) |
| `size_t region` | `1` | `1`: point-wise QoI. `r > 1`: regional QoI over blocks of `r` values per dimension. See [Concepts](concepts.md). |
| `double qoi_tol` | `1e-3` | QoI error threshold |
| `BoundMode qoi_tol_mode` | `REL` | `BoundMode::ABS` or `BoundMode::REL` (relative to the value range of the QoI) |
| `double data_eb` | `0` | bound on the error of the data; `0` means none |
| `BoundMode data_eb_mode` | `REL` | `ABS` or `REL` (relative to the value range of the data) |
| `int threads` | `0` | threads for QPET's own passes and for backends whose output does not depend on it; `0` means all |
| `bool backend_parallel` | `false` | let SZ3 and SZo use their OpenMP mode, which is faster and compresses less |
| `double quantile` | `-1` | if `>= 0`, skip the tuning and use this quantile of the estimated point-wise error bounds |
| `double backend_eb` | `0` | if `> 0`, skip the tuning and run the backend with this absolute error bound |
| `TuningParams tuning` | all unset | how the tuning samples and searches, see [Tuning](tuning.md) |
| `std::map<std::string, std::string> backend_options` | empty | options for the backend, see [Backends](backends.md) |
| `bool verbose` | `false` | print what the tuning tried to standard error |

`conf.set(name, value)` sets a field from text, which is convenient when options come from a file or a command line:
`set("qoi_tol_mode", "abs")`, `set("region", "4")`, `set("tuning.min_blocks", "24")`, `set("backend.algo", "interp")`.
It throws `std::invalid_argument` for an unknown name or a bad value.

### `CompressInfo`

| field | meaning |
|---|---|
| `qoi_tol_abs`, `data_eb_abs` | the absolute tolerances in effect (`data_eb_abs` is `+inf` without a data bound) |
| `backend_eb` | the error bound the backend ran with |
| `quantile` | the share of the sampled values whose estimated error bound is below `backend_eb` |
| `candidates` | error bounds the tuning tried |
| `num` | number of values |
| `outliers` | values stored losslessly |
| `base_bytes`, `patch_bytes`, `total_bytes` | size of the backend's stream, of the outliers, of the whole stream |
| `max_qoi_err`, `max_data_err` | the largest errors after correction, as the compression measured them |
| `t_estimate`, `t_tune`, `t_compress`, `t_decompress`, `t_validate`, `t_total` | seconds per stage |

## Inspect and decompress

```cpp
StreamInfo inspect(const uint8_t *buf, size_t size);
void decompress(const uint8_t *buf, size_t size, void *out, int threads = 0);
```

`inspect` reads the header. `decompress` writes `StreamInfo::num` values of `StreamInfo::dtype` to `out`, which the
caller allocates. The backend named in the stream must be part of the build.

`StreamInfo` has `backend`, `backend_version`, `qoi`, `dtype`, `dims`, `num`, `region`, `qoi_tol_abs`, `data_eb_abs`,
`backend_eb` and `outliers`.

## Verify

```cpp
VerifyResult verify(const void *orig, const void *dec, DType dtype, const std::vector<size_t> &dims,
                    const std::string &qoi, size_t region = 1, int threads = 0);
VerifyResult verify(const void *orig, const void *dec, DType dtype, size_t num, const std::string &qoi, int threads = 0);
```

Computes the errors between two arrays, without any knowledge of how one was made from the other. The second form is
a point-wise QoI on a flat array.

| field | meaning |
|---|---|
| `max_qoi_err` | largest QoI error: point-wise, or of the regional means with `region > 1`. Infinite if a value where the QoI is not finite was not reproduced exactly. |
| `max_data_err` | largest `\|x - x'\|` |
| `qoi_range`, `data_range` | value ranges over the original, the ones a relative tolerance refers to |
| `psnr` | of the data, in dB; infinite if the arrays are equal |

## Example

```cpp
#include <vector>
#include "qpet/qpet.hpp"

std::vector<float> data = ...;                       // nz * ny * nx values
qpet::Config conf;
conf.backend = "sperr";
conf.qoi = "x^2";
conf.region = 4;                                     // mean of x^2 over 4x4x4 blocks
conf.qoi_tol = 1e-3;                                 // relative to the range of those means
conf.data_eb = 1e-2;                                 // and every value within 1 % of the value range

qpet::CompressInfo info;
std::vector<uint8_t> cmp = qpet::compress(data.data(), {nz, ny, nx}, conf, &info);

qpet::StreamInfo s = qpet::inspect(cmp.data(), cmp.size());
std::vector<float> dec(s.num);
qpet::decompress(cmp.data(), cmp.size(), dec.data());

qpet::VerifyResult v = qpet::verify(data.data(), dec.data(), s.dtype, s.dims, s.qoi, s.region);
// v.max_qoi_err <= info.qoi_tol_abs and v.max_data_err <= info.data_eb_abs
```

## Threads

`compress`, `decompress` and `verify` use OpenMP inside and can be called from several threads at once on different
data. `Config::threads` limits the threads of one call.

## Lower-level pieces

- `qpet::Expr` (`qpet/expr.hpp`): a compiled QoI expression. `eval(x)`, `eval(const double *x, double *out, size_t n)`
  for arrays, `jet(x)` for the value with the first two derivatives, `is_linear()`.
- `qpet::pointwise_eb(f, x, tau, cap)`: the estimated error bound of one value.
- `qpet::Backend`, `get_backend(name)`, `available_backends()` (`qpet/backend.hpp`): the compressors, see
  [Backends](backends.md).
