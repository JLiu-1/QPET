# QPET documentation

| page | what it covers |
|---|---|
| [Concepts](concepts.md) | what QPET guarantees, point-wise and regional QoIs, tolerances, how the compression works, limits |
| [Building](build.md) | requirements, CMake options, how the compressors are obtained, what the build produces, tests, use from another project |
| [Command-line tool](cli.md) | every option of `qpet`, the report it prints, exit codes, examples |
| [QoI expressions](qoi.md) | the syntax of the QoI function, how it is evaluated, what to watch for |
| [C++ API](cpp-api.md) | `qpet::compress`, `decompress`, `inspect`, `verify`, the configuration and result structures, `Expr` |
| [C API](c-api.md) | the functions of `libqpet_c` |
| [Python API](python.md) | the `qpet` package |
| [Tuning](tuning.md) | how the backend's error bound is chosen, every tuning parameter, how to diagnose a choice |
| [Backends](backends.md) | SZ3, SZo, ZFP and SPERR: modes, limits, options; how to add your own compressor |
| [Stream format](format.md) | the layout of a compressed stream |

If you only want to try QPET, the [README](../README.md) has the build command and a first example.
