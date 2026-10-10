# Python API

The package `qpet` works on NumPy arrays. It is a thin layer over the [C API](c-api.md): it loads the shared library
`libqpet_c` that the build produces.

## Getting it

Build QPET (see [Building](build.md)). The build puts the package, with the library next to it, in `build/python`.
Either point Python at it or install it:

```bash
export PYTHONPATH=/path/to/QPET/build/python      # use it in place
# or
pip install /path/to/QPET/build/python            # install it (copies the library into the environment)
```

It needs Python 3.8 or newer and NumPy. The library is looked for, in this order: the path in the environment variable
`QPET_LIBRARY`, next to the package, the system's library path.

The package is tied to the build it came from: it has the backends of that build, and the library was compiled for that
machine. After rebuilding QPET, an installed copy has to be installed again.

## Functions

Arrays are `float32` or `float64` with 1 to 4 dimensions. Anything NumPy can view as such is accepted; a
non-contiguous array or view is copied first. The shape is used as it is (C order).

```python
qpet.compress(data, qoi, qoi_tol, *, qoi_tol_mode="rel", region=1, data_eb=0.0, data_eb_mode="rel",
              backend="sz3", threads=0, backend_parallel=False, quantile=None, backend_eb=None,
              tuning=None, backend_options=None, verbose=False, return_info=False)
```

Returns the compressed stream as `bytes`, or `(bytes, dict)` with `return_info=True`.

| argument | meaning |
|---|---|
| `data` | the array |
| `qoi` | the QoI function as text, see [QoI expressions](qoi.md) |
| `qoi_tol`, `qoi_tol_mode` | QoI error threshold; `"rel"` (relative to the value range of the QoI) or `"abs"` |
| `region` | `1`: point-wise QoI. `r > 1`: regional QoI over blocks of `r` values per dimension. See [Concepts](concepts.md). |
| `data_eb`, `data_eb_mode` | bound on the error of the data, `0` for none; `"rel"` or `"abs"` |
| `backend` | `"sz3"`, `"szo"`, `"zfp"` or `"sperr"` |
| `threads` | `0` for all |
| `backend_parallel` | let SZ3 and SZo use their OpenMP mode (faster, lower compression ratio) |
| `quantile`, `backend_eb` | skip the tuning, see [Tuning](tuning.md) |
| `tuning` | dict of [tuning parameters](tuning.md), e.g. `{"min_blocks": 24, "quantiles": [0, 0.1, 1]}` |
| `backend_options` | dict of options for the [backend](backends.md), e.g. `{"algo": "interp"}` |
| `verbose` | print what the tuning tried to standard error |

The info dict has `compression_ratio` and the fields listed for `CompressInfo` in the [C++ API](cpp-api.md):
`qoi_tol_abs`, `data_eb_abs`, `backend_eb`, `quantile`, `candidates`, `num`, `outliers`, `base_bytes`, `patch_bytes`,
`total_bytes`, `max_qoi_err`, `max_data_err` and the timings.

```python
qpet.decompress(blob, threads=0)
```

Returns an array with the shape and dtype of the original. `blob` can be `bytes`, a `bytearray`, a `memoryview` or a
`uint8` array.

```python
qpet.inspect(blob)
```

Returns a dict: `backend`, `backend_version`, `qoi`, `region`, `shape`, `dtype`, `num`, `qoi_tol_abs`, `data_eb_abs`,
`backend_eb`, `outliers`.

```python
qpet.verify(original, decompressed, qoi, region=1, threads=0)
```

Returns a dict: `max_qoi_err`, `max_data_err`, `qoi_range`, `data_range`, `psnr`, computed from the two arrays alone.

```python
qpet.backends()
```

Returns the backends of the build and the versions of their compressors, e.g. `{"sz3": "3.4.0", "zfp": "1.0.1"}`.

## Errors

A failure reported by the library raises `qpet.QpetError` (a `RuntimeError`) with the library's message: an unknown
backend, a QoI that does not parse, an unknown tuning parameter, a stream that is not a QPET stream. An array of the
wrong type raises `TypeError`, one with too many dimensions `ValueError`.

## Example

```python
import numpy as np
import qpet

data = np.fromfile("Uf48.bin.f32", dtype=np.float32).reshape(100, 500, 500)

blob, info = qpet.compress(data, qoi="x^2", qoi_tol=1e-3, data_eb=1e-2, backend="sperr", return_info=True)
print(info["compression_ratio"], info["outliers"])
open("U.qpet", "wb").write(blob)

back = qpet.decompress(blob)                       # shape (100, 500, 500), float32
errors = qpet.verify(data, back, qoi="x^2")
assert errors["max_qoi_err"] <= info["qoi_tol_abs"]

# the mean of x^2 over 4x4x4 blocks instead
blob = qpet.compress(data, qoi="x^2", qoi_tol=1e-3, region=4, backend="sz3")
print(qpet.inspect(blob))
```

The functions release the interpreter lock only as far as `ctypes` does for any foreign call, which is for the whole
duration of the call: several Python threads can compress different arrays at the same time.
