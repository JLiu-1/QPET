"""Round trips through the Python package. Run with PYTHONPATH=<build>/python."""
import numpy as np

import qpet

failures = 0


def check(cond, what):
    global failures
    if not cond:
        failures += 1
        print("FAIL", what)


rng = np.random.default_rng(7)
z, y, x = np.meshgrid(np.arange(40), np.arange(56), np.arange(72), indexing="ij")
field = (2 + np.sin(0.07 * z) / 3 + np.sin(0.14 * y + 1) / 3 + np.sin(0.21 * x + 2) / 3
         + rng.uniform(-0.005, 0.005, z.shape))

names = qpet.backends()
check(len(names) > 0, "no backends")
for backend in names:
    for dtype in (np.float32, np.float64):
        for region in (1, 4):
            data = field.astype(dtype)
            blob, info = qpet.compress(data, qoi="x^2", qoi_tol=1e-3, data_eb=1e-2, backend=backend, region=region,
                                       threads=2, return_info=True)
            back = qpet.decompress(blob)
            head = qpet.inspect(blob)
            v = qpet.verify(data, back, qoi="x^2", region=region)
            what = f"{backend} {np.dtype(dtype)} region={region}"
            check(isinstance(blob, bytes) and len(blob) == info["total_bytes"], what + ": stream")
            check(back.shape == data.shape and back.dtype == data.dtype, what + ": shape and dtype")
            check(head["backend"] == backend and head["qoi"] == "x^2" and head["region"] == region
                  and head["shape"] == data.shape and head["dtype"] == data.dtype, what + ": inspect")
            check(v["max_qoi_err"] <= info["qoi_tol_abs"] and v["max_data_err"] <= info["data_eb_abs"], what + ": bounds")
            check(info["compression_ratio"] > 1, what + ": ratio")

# a non-contiguous view, absolute tolerance, tuning parameters and a backend option
first = next(iter(names))
view = field.astype(np.float32)[::2, :, ::3]
blob = qpet.compress(view, qoi="log2(x)", qoi_tol=1e-3, qoi_tol_mode="abs", backend=first,
                     tuning={"sample_rate": 0.05, "quantiles": [0, 0.1, 1]})
check(qpet.verify(view, qpet.decompress(blob), qoi="log2(x)")["max_qoi_err"] <= 1e-3, "absolute tolerance on a view")
check(np.array_equal(qpet.decompress(bytearray(blob)), qpet.decompress(memoryview(blob))), "bytes-like streams")

for bad in (lambda: qpet.compress(field.astype(np.float32), qoi="nonsense(x)", qoi_tol=1e-3),
            lambda: qpet.compress(field.astype(np.float32), qoi="x", qoi_tol=1e-3, backend="nope"),
            lambda: qpet.compress(field.astype(np.float32), qoi="x", qoi_tol=1e-3, tuning={"nope": 1}),
            lambda: qpet.decompress(b"not a stream")):
    try:
        bad()
        check(False, "an error was not reported")
    except qpet.QpetError:
        pass
try:
    qpet.compress(field.astype(np.int32), qoi="x", qoi_tol=1e-3)
    check(False, "an integer array was accepted")
except TypeError:
    pass

print(f"{failures} FAILURES" if failures else "all Python tests passed")
raise SystemExit(1 if failures else 0)
