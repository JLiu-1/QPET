# C API

```c
#include "qpet/qpet_c.h"
```

The C interface is the shared library `libqpet_c`, which contains QPET and the compressors, so nothing else has to be
linked. In CMake, link the target `qpet_c`.

Conventions:

- Every function returns `0` on success. On failure it returns non-zero and `qpet_last_error()` gives the reason. The
  text belongs to the calling thread and is valid until its next call into the library.
- Dimensions are in C order (slowest-varying first), 1 to 4 of them.
- `dtype` is `QPET_F32` or `QPET_F64`.
- Buffers and texts that a function hands out are released with `qpet_free()`.
- Text results are lines of `name=value`.

## Functions

```c
int qpet_compress(const void *data, int dtype, size_t ndims, const size_t *dims, const char *options,
                  unsigned char **out, size_t *out_size, char **info);
```

Compresses. `options` is a text with one `name=value` per line; `NULL` or `""` keeps every default. The names are the
fields of the C++ `Config` (see [C++ API](cpp-api.md)):

| name | value |
|---|---|
| `backend` | `sz3`, `szo`, `zfp`, `sperr` |
| `qoi` | the QoI function |
| `region` | region edge, `1` for a point-wise QoI |
| `qoi_tol`, `data_eb` | numbers |
| `qoi_tol_mode`, `data_eb_mode` | `abs` or `rel` |
| `threads` | integer, `0` for all |
| `backend_parallel`, `verbose` | `0` or `1` |
| `quantile`, `backend_eb` | numbers; either one skips the tuning |
| `tuning.<parameter>` | a [tuning parameter](tuning.md), e.g. `tuning.min_blocks=24` |
| `backend.<option>` | an option of the [backend](backends.md), e.g. `backend.algo=interp` |

`*out` and `*out_size` receive the stream. If `info` is not `NULL`, `*info` receives the fields of the C++
`CompressInfo`: `qoi_tol_abs`, `data_eb_abs`, `backend_eb`, `quantile`, `candidates`, `num`, `outliers`, `base_bytes`,
`patch_bytes`, `total_bytes`, `max_qoi_err`, `max_data_err`, and the timings `t_estimate`, `t_tune`, `t_compress`,
`t_decompress`, `t_validate`, `t_total`.

```c
int qpet_inspect(const unsigned char *buf, size_t size, char **info);
```

Reads the header of a stream: `backend`, `backend_version`, `qoi`, `region`, `dtype` (`0` or `1`), `dims`
(comma-separated), `num`, `qoi_tol_abs`, `data_eb_abs`, `backend_eb`, `outliers`.

```c
int qpet_decompress(const unsigned char *buf, size_t size, void *out, size_t out_bytes, int threads);
```

Decompresses into `out`, which the caller allocates: `num` values of the stream's type. `out_bytes` is its size in
bytes and must match the stream, otherwise the call fails without writing.

```c
int qpet_verify(const void *orig, const void *dec, int dtype, size_t ndims, const size_t *dims, const char *qoi,
                size_t region, int threads, char **result);
```

Errors between two arrays: `max_qoi_err`, `max_data_err`, `qoi_range`, `data_range`, `psnr`.

```c
int qpet_backends(char **list);       /* lines of name=version */
const char *qpet_last_error(void);
void qpet_free(void *p);
```

## Example

```c
#include <stdio.h>
#include <stdlib.h>
#include "qpet/qpet_c.h"

int roundtrip(const float *data, size_t nz, size_t ny, size_t nx) {
    size_t dims[3] = {nz, ny, nx}, size = 0;
    unsigned char *cmp = NULL;
    char *info = NULL;
    const char *options = "backend=sz3\nqoi=x^2\nqoi_tol=1e-3\ndata_eb=1e-2\n";
    if (qpet_compress(data, QPET_F32, 3, dims, options, &cmp, &size, &info)) {
        fprintf(stderr, "%s\n", qpet_last_error());
        return 1;
    }
    printf("%s", info);                                  /* name=value lines */

    float *back = malloc(nz * ny * nx * sizeof(float));
    int rc = qpet_decompress(cmp, size, back, nz * ny * nx * sizeof(float), 0);
    if (rc) fprintf(stderr, "%s\n", qpet_last_error());

    free(back);
    qpet_free(info);
    qpet_free(cmp);
    return rc;
}
```
