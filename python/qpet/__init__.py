"""QPET: QoI-preserving error-bounded lossy compression on top of SZ3, SZo, ZFP or SPERR.

    import numpy as np, qpet

    data = np.fromfile("density.f32", dtype=np.float32).reshape(256, 384, 384)
    blob = qpet.compress(data, qoi="x^2", qoi_tol=1e-3, backend="sz3")
    back = qpet.decompress(blob)          # same shape and dtype as data
    print(qpet.verify(data, back, qoi="x^2"))

Arrays are float32 or float64 with 1 to 4 dimensions, in C order.
"""
import ctypes
import ctypes.util
import os
import sys

import numpy as np

__all__ = ["compress", "decompress", "inspect", "verify", "backends", "QpetError"]
__version__ = "2.0.0"


class QpetError(RuntimeError):
    """A failure reported by the QPET library."""


def _load():
    names = {"darwin": ["libqpet_c.dylib"], "win32": ["qpet_c.dll", "libqpet_c.dll"]}.get(sys.platform, ["libqpet_c.so"])
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [os.environ["QPET_LIBRARY"]] if os.environ.get("QPET_LIBRARY") else []
    candidates += [os.path.join(here, n) for n in names]
    found = ctypes.util.find_library("qpet_c")
    if found:
        candidates.append(found)
    for path in candidates:
        if os.environ.get("QPET_LIBRARY") == path or os.path.exists(path) or path == found:
            try:
                return ctypes.CDLL(path)
            except OSError as e:
                raise ImportError(f"cannot load the QPET library {path}: {e}") from None
    raise ImportError("cannot find the QPET library (libqpet_c). Build QPET and use the package under <build>/python, "
                      "or set QPET_LIBRARY to the path of the library.")


_lib = _load()
_ubyte_p = ctypes.POINTER(ctypes.c_ubyte)
_size_p = ctypes.POINTER(ctypes.c_size_t)
_text_p = ctypes.POINTER(ctypes.c_void_p)
_lib.qpet_compress.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_size_t, _size_p, ctypes.c_char_p,
                               ctypes.POINTER(_ubyte_p), _size_p, _text_p]
_lib.qpet_inspect.argtypes = [ctypes.c_void_p, ctypes.c_size_t, _text_p]
_lib.qpet_decompress.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
_lib.qpet_verify.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_size_t, _size_p,
                             ctypes.c_char_p, ctypes.c_size_t, ctypes.c_int, _text_p]
_lib.qpet_backends.argtypes = [_text_p]
_lib.qpet_last_error.restype = ctypes.c_char_p
_lib.qpet_free.argtypes = [ctypes.c_void_p]

_DTYPES = {np.dtype(np.float32): 0, np.dtype(np.float64): 1}


def _check(rc):
    if rc != 0:
        raise QpetError(_lib.qpet_last_error().decode())


def _take_text(ptr):
    """Reads and frees a "name=value" text the library handed out."""
    try:
        text = ctypes.string_at(ptr).decode()
    finally:
        _lib.qpet_free(ptr)
    out = {}
    for line in text.splitlines():
        name, _, value = line.partition("=")
        out[name] = value
    return out


def _numbers(d, keep=()):
    for k, v in d.items():
        if k in keep:
            continue
        try:
            d[k] = int(v)
        except ValueError:
            d[k] = float(v)
    return d


def _array(data, what="data"):
    a = np.ascontiguousarray(data)
    if a.dtype not in _DTYPES:
        raise TypeError(f"{what} must be float32 or float64, not {a.dtype}")
    if not 1 <= a.ndim <= 4:
        raise ValueError(f"{what} must have 1 to 4 dimensions, not {a.ndim}")
    return a


def _dims(a):
    return (ctypes.c_size_t * a.ndim)(*a.shape)


def _stream(blob):
    return np.frombuffer(blob, dtype=np.uint8)


def compress(data, qoi, qoi_tol, *, qoi_tol_mode="rel", region=1, data_eb=0.0, data_eb_mode="rel", backend="sz3",
             threads=0, backend_parallel=False, quantile=None, backend_eb=None, tuning=None, backend_options=None,
             verbose=False, return_info=False):
    """Compresses an array so that the QoI computed from the decompressed array stays within a tolerance.

    data             float32 or float64 array with 1 to 4 dimensions
    qoi              point-wise QoI f(x) as text, e.g. "x^2", "log2(x)", "sin(10x)"
    qoi_tol          QoI error threshold
    qoi_tol_mode     "rel": relative to the value range of the QoI over the data; "abs": absolute
    region           1: bound the error of f at every value. r > 1: bound the error of the mean of f over every block
                     of r values per dimension
    data_eb          bound on the error of the data itself, 0 for none; data_eb_mode as qoi_tol_mode
    backend          "sz3", "szo", "zfp" or "sperr" (see backends())
    threads          0: all
    backend_parallel let SZ3/SZo use their OpenMP mode (faster, lower compression ratio)
    quantile         skip the tuning: use this quantile of the point-wise error bounds as the backend's error bound
    backend_eb       skip the tuning: run the backend with this absolute error bound
    tuning           dict of tuning parameters, e.g. {"sample_rate": 0.02, "min_blocks": 24}
    backend_options  dict of options for the backend, e.g. {"algo": "interp"}
    verbose          print what the tuning tried to stderr
    return_info      also return a dict with what the compression did (sizes, outliers, errors, timings)

    Returns the compressed stream as bytes, or (bytes, dict) with return_info.
    """
    a = _array(data)
    options = {"backend": backend, "qoi": qoi, "qoi_tol": repr(float(qoi_tol)), "qoi_tol_mode": qoi_tol_mode,
               "region": int(region), "data_eb": repr(float(data_eb)), "data_eb_mode": data_eb_mode,
               "threads": int(threads), "backend_parallel": int(bool(backend_parallel)), "verbose": int(bool(verbose))}
    if quantile is not None:
        options["quantile"] = repr(float(quantile))
    if backend_eb is not None:
        options["backend_eb"] = repr(float(backend_eb))
    for k, v in (tuning or {}).items():
        if k == "quantiles" and not isinstance(v, str):
            v = ":".join(repr(float(q)) for q in v)
        elif isinstance(v, bool):
            v = int(v)
        options["tuning." + k] = v
    for k, v in (backend_options or {}).items():
        options["backend." + k] = v
    text = "".join(f"{k}={v}\n" for k, v in options.items()).encode()

    out, size, info = _ubyte_p(), ctypes.c_size_t(), ctypes.c_void_p()
    _check(_lib.qpet_compress(a.ctypes.data, _DTYPES[a.dtype], a.ndim, _dims(a), text, ctypes.byref(out),
                              ctypes.byref(size), ctypes.byref(info)))
    try:
        blob = ctypes.string_at(out, size.value)
    finally:
        _lib.qpet_free(out)
    stats = _numbers(_take_text(info))
    if return_info:
        stats["compression_ratio"] = a.nbytes / len(blob)
        return blob, stats
    return blob


def inspect(blob):
    """What a compressed stream holds: backend, QoI, region, dtype, shape, tolerances, number of outliers."""
    s = _stream(blob)
    info = ctypes.c_void_p()
    _check(_lib.qpet_inspect(s.ctypes.data, s.size, ctypes.byref(info)))
    d = _take_text(info)
    shape = tuple(int(x) for x in d.pop("dims").split(","))
    dtype = np.dtype(np.float32 if d.pop("dtype") == "0" else np.float64)
    d = _numbers(d, keep=("backend", "backend_version", "qoi"))
    d["shape"], d["dtype"] = shape, dtype
    return d


def decompress(blob, threads=0):
    """Decompresses a stream made by compress(). The array has the shape and dtype of the original."""
    info = inspect(blob)
    s = _stream(blob)
    out = np.empty(info["shape"], dtype=info["dtype"])
    _check(_lib.qpet_decompress(s.ctypes.data, s.size, out.ctypes.data, out.nbytes, int(threads)))
    return out


def verify(original, decompressed, qoi, region=1, threads=0):
    """Errors between two arrays: max_qoi_err, max_data_err, qoi_range, data_range, psnr.

    With region > 1, max_qoi_err and qoi_range are those of the regional QoI (the mean of f over each block).
    """
    a, b = _array(original, "original"), _array(decompressed, "decompressed")
    if a.shape != b.shape or a.dtype != b.dtype:
        raise ValueError("the two arrays must have the same shape and dtype")
    result = ctypes.c_void_p()
    _check(_lib.qpet_verify(a.ctypes.data, b.ctypes.data, _DTYPES[a.dtype], a.ndim, _dims(a), qoi.encode(),
                            int(region), int(threads), ctypes.byref(result)))
    return _numbers(_take_text(result))


def backends():
    """The backends of this build and the versions of their compressors, e.g. {"sz3": "3.4.0", ...}."""
    text = ctypes.c_void_p()
    _check(_lib.qpet_backends(ctypes.byref(text)))
    return _take_text(text)
