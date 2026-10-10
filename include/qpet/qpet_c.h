/*
 * C interface of QPET, in the shared library qpet_c. It is also what the Python package calls.
 *
 * Every function returns 0 on success. On failure it returns a non-zero value and qpet_last_error() gives the reason.
 * Dimensions are in C order (slowest-varying first). Text results are lines of "name=value"; buffers and texts that a
 * function hands out are released with qpet_free().
 */
#ifndef QPET_C_H
#define QPET_C_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define QPET_F32 0
#define QPET_F64 1

/*
 * options: lines of "name=value" (the names of qpet::Config):
 *     backend, qoi, region, qoi_tol, qoi_tol_mode (abs | rel), data_eb, data_eb_mode (abs | rel), threads,
 *     backend_parallel (0 | 1), quantile, backend_eb, verbose (0 | 1),
 *     tuning.<parameter> (see qpet::TuningParams), backend.<option> (passed to the backend).
 *   NULL or "" keeps every default.
 * out, out_size: the compressed stream.
 * info: if not NULL, receives what the compression did (the fields of qpet::CompressInfo).
 */
int qpet_compress(const void *data, int dtype, size_t ndims, const size_t *dims, const char *options,
                  unsigned char **out, size_t *out_size, char **info);

/* info: backend, backend_version, qoi, region, dtype (0 | 1), dims (comma-separated), num, the tolerances, outliers */
int qpet_inspect(const unsigned char *buf, size_t size, char **info);

/* out holds num values of the stream's dtype; out_bytes is its size and is checked against the stream */
int qpet_decompress(const unsigned char *buf, size_t size, void *out, size_t out_bytes, int threads);

/* result: max_qoi_err, max_data_err, qoi_range, data_range, psnr */
int qpet_verify(const void *orig, const void *dec, int dtype, size_t ndims, const size_t *dims, const char *qoi,
                size_t region, int threads, char **result);

/* lines of "name=version", one per backend of this build */
int qpet_backends(char **list);

/* reason of the last failure in this thread; valid until the next call into the library from this thread */
const char *qpet_last_error(void);

void qpet_free(void *p);

#ifdef __cplusplus
}
#endif

#endif
