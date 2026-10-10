# Stream format

A QPET stream is self-describing: decompression needs nothing else. All integers and floating-point numbers are
little-endian; QPET does not convert byte order, so a stream is read on a machine of the same endianness as the one
that wrote it.

```
"QPET"                      4 bytes
format                      u8        2 (1 is the same without the region field, and is still read)
dtype                       u8        0: float32, 1: float64
ndims                       u8        1 to 4
reserved                    u8        0
dims[ndims]                 u64 each  C order, slowest-varying first
region                      u64       1: point-wise QoI, otherwise the region edge       (format 2)
backend_eb                  f64       absolute error bound the backend ran with
qoi_tol_abs                 f64       absolute QoI tolerance
data_eb_abs                 f64       absolute data bound, +inf if none
backend                     string    name of the backend
backend_version             string    version of the compressor that wrote the backend's stream
qoi                         string    the QoI expression
base_bytes                  u64
backend stream              base_bytes bytes, as the compressor produced it
outliers                    u64       number of values stored losslessly
  if outliers > 0:
  index gaps                zstd block
  values                    zstd block
```

- **string**: `u16` length, then the bytes, without a terminator.
- **zstd block**: `u64` uncompressed size, `u64` compressed size, then a Zstandard frame.
- **index gaps**: the flat indices of the outliers in increasing order, each stored as the difference from the
  previous one (the first from 0), as an unsigned LEB128 integer.
- **values**: the original values of the outliers, in the same order, in the stream's data type.

To decompress: decompress the backend stream with the named backend (`backend_eb` is passed along, because ZFP's
fixed-accuracy mode needs the tolerance to decode), then write each outlier value at its index.

The QoI, the tolerances and the region are recorded for reference and are not needed to decompress.

## Compatibility

- The backend's stream is whatever the compressor produced. Whether a later version of that compressor can read it
  is the compressor's own compatibility policy; `backend_version` records what wrote it.
- Decompression fails with a clear error if the stream names a backend that the build does not have.
- A stream is rejected if its magic, format number or sizes are inconsistent. There is no checksum.
