# Command-line tool

`qpet` compresses a raw binary file of `float32` or `float64` values, and decompresses a QPET stream. Run it without
arguments for a summary of the options.

## Compress

```
qpet -b <backend> -f|-d -<N> <n1> ... <nN> -i <input> -q <QoI> -m <mode> <tol>
     [-r <edge>] [-M <mode> <bound>] [-z <compressed>] [-o <decompressed>] [-a] [options]
```

| option | meaning |
|---|---|
| `-b <backend>` | `sz3`, `szo`, `zfp` or `sperr` (default `sz3`). `qpet -l` lists the backends of the build with their versions. See [Backends](backends.md). |
| `-f`, `-d` | the input is `float32`, `float64`. Required. |
| `-1 nx`, `-2 nx ny`, `-3 nx ny nz`, `-4 nx ny nz nw` | dimensions, **fastest-varying first**: a C array `a[100][200][300]` is `-3 300 200 100`. Required. This is the convention of the SZ3 and SPERR tools. |
| `-i <file>` | input file. Its size must be exactly the product of the dimensions times 4 or 8 bytes. |
| `-q <QoI>` | the QoI function, see [QoI expressions](qoi.md). Quote it. Default `x`. |
| `-r <edge>` | regional QoI: bound the error of the mean of the QoI over blocks of `<edge>` values per dimension, instead of the error at every value. See [Concepts](concepts.md). |
| `-m ABS\|REL <tol>` | QoI error threshold. Required. `REL` is relative to the value range of the QoI. |
| `-M ABS\|REL <bound>` | bound on the error of the data itself. `REL` is relative to the value range of the data. Default: none. |
| `-z <file>` | write the compressed stream. Without it the stream is only measured. |
| `-o <file>` | decompress the stream and write the result. |
| `-a` | decompress and report the errors, computed from the original and the decompressed data independently of the compression. |
| `-t <n>` | threads. Default: all. |
| `-P` | let SZ3 and SZo run in their OpenMP mode. They are faster and compress less in that mode. |
| `-O <key>=<value>` | option for the backend, for example `-O algo=interp` (`sz3`, `szo`) or `-O chunk=128` (`sperr`). Repeatable. See [Backends](backends.md). |
| `-s <rate>` | fraction of the input the tuning compresses per candidate (default 0.01). Same as `-S sample_rate=<rate>`. |
| `-S <name>=<value>` | tuning parameter, see [Tuning](tuning.md). Repeatable. |
| `-Q <quantile>` | skip the tuning: use this quantile (0 to 1) of the estimated point-wise error bounds as the backend's error bound. |
| `-E <bound>` | skip the tuning: run the backend with this absolute error bound. |
| `-v` | print every candidate the tuning tried, to standard error. |

## Decompress

```
qpet -z <compressed> -o <decompressed> [-t <n>]
```

The stream records the backend, the data type and the dimensions, so none of them is given. The backend must be part
of the build.

## Other

| option | meaning |
|---|---|
| `-l` | list the backends of this build and the versions of their compressors |
| `-h`, `--help` | usage |

## The report

Compression prints one `name = value` line per item:

| line | meaning |
|---|---|
| `backend` | backend and the version of its compressor |
| `QoI`, `QoI region` | the QoI function, and the region edge if the QoI is regional |
| `QoI tolerance (abs)` | the absolute QoI tolerance in effect |
| `data error bound (abs)` | the absolute data bound in effect, `inf` if none |
| `backend error bound` | the error bound the backend ran with |
| `backend error bound quantile` | the share of the sampled values whose estimated error bound is below it |
| `tuning candidates` | how many error bounds the tuning tried (0 if it was skipped or not needed) |
| `outliers` | values stored losslessly, and their share |
| `backend bytes`, `outlier bytes`, `compressed bytes` | size of the backend's stream, of the outliers, of the whole stream |
| `compression ratio`, `bit rate` | input bytes / compressed bytes; compressed bits per value |
| `time estimate`, `time tune`, `time backend compress`, `time backend decompress`, `time validate` | seconds per stage |
| `compression time`, `compression speed` | total seconds, and input megabytes (10^6 bytes) per second. File reading and writing are not included. |

With `-a` or `-o`:

| line | meaning |
|---|---|
| `decompression time`, `decompression speed` | |
| `max data error` | largest `\|x - x'\|`, absolute and relative to the value range (`-a`) |
| `max QoI error` | largest QoI error (point-wise, or of the regional means), absolute and relative to the range of the QoI (`-a`) |
| `PSNR` | of the data, in dB (`-a`) |
| `bounds respected` | `yes` if both maxima are within the tolerances in effect (`-a`) |

## Exit codes

| code | meaning |
|---|---|
| 0 | success |
| 1 | error (bad arguments, unreadable file, failure in a backend); the reason is on standard error |
| 2 | `-a` found a bound broken. This is not expected to happen; please report it. |

## Examples

```bash
# point-wise x^2 within 0.1 % of its range, data within 1 % of its range, with SZ3; report the errors
qpet -b sz3 -f -3 500 500 100 -i Uf48.bin.f32 -q "x^2" -m REL 1e-3 -M REL 1e-2 -z U.qpet -a

# the mean of x^2 over 4x4x4 blocks instead, with SPERR
qpet -b sperr -f -3 500 500 100 -i Uf48.bin.f32 -q "x^2" -r 4 -m REL 1e-3 -M REL 1e-2 -z U.qpet -a

# an absolute tolerance on log2, double precision, 2-D
qpet -b zfp -d -2 3600 1800 -i field.f64 -q "log2(x)" -m ABS 0.01 -z field.qpet

# decompress
qpet -z U.qpet -o U.out.f32

# see what the tuning does, with more sample blocks
qpet -b sz3 -f -2 3600 1800 -i field.f32 -q "x^3" -m REL 1e-2 -S min_blocks=48 -v
```
