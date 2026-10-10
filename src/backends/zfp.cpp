#include <stdexcept>

#include "qpet/backend.hpp"
#include "zfp.h"

namespace qpet {

namespace {

zfp_field *make_field(void *data, DType dtype, const std::vector<size_t> &dims) {
    zfp_type type = dtype == DType::F32 ? zfp_type_float : zfp_type_double;
    const size_t n = dims.size();  // zfp takes the fastest-varying dimension first
    switch (n) {
        case 1: return zfp_field_1d(data, type, dims[0]);
        case 2: return zfp_field_2d(data, type, dims[1], dims[0]);
        case 3: return zfp_field_3d(data, type, dims[2], dims[1], dims[0]);
        default: return zfp_field_4d(data, type, dims[3], dims[2], dims[1], dims[0]);
    }
}

/** fixed-accuracy mode */
class ZFPBackend : public Backend {
   public:
    std::string name() const override { return "zfp"; }

    std::string version() const override { return ZFP_VERSION_STRING; }

    bool supports(DType, size_t ndims) const override { return ndims >= 1 && ndims <= 4; }

    std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims,
                                  const BackendParams &params) const override {
        zfp_field *field = make_field(const_cast<void *>(data), dtype, dims);
        zfp_stream *zfp = zfp_stream_open(nullptr);
        zfp_stream_set_accuracy(zfp, params.abs_eb);
        if (params.threads > 1 && zfp_stream_set_execution(zfp, zfp_exec_omp))
            zfp_stream_set_omp_threads(zfp, static_cast<uint>(params.threads));
        std::vector<uint8_t> out(zfp_stream_maximum_size(zfp, field));
        bitstream *bits = stream_open(out.data(), out.size());
        zfp_stream_set_bit_stream(zfp, bits);
        zfp_stream_rewind(zfp);
        size_t size = zfp_compress(zfp, field);
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(bits);
        if (size == 0) throw std::runtime_error("QPET: zfp compression failed");
        out.resize(size);
        return out;
    }

    void decompress(const uint8_t *buf, size_t size, DType dtype, const std::vector<size_t> &dims,
                    const BackendParams &params, void *out) const override {
        zfp_field *field = make_field(out, dtype, dims);
        zfp_stream *zfp = zfp_stream_open(nullptr);
        zfp_stream_set_accuracy(zfp, params.abs_eb);
        bitstream *bits = stream_open(const_cast<uint8_t *>(buf), size);
        zfp_stream_set_bit_stream(zfp, bits);
        zfp_stream_rewind(zfp);
        size_t read = zfp_decompress(zfp, field);
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(bits);
        if (read == 0) throw std::runtime_error("QPET: zfp decompression failed");
    }

    /**
     * zfp codes every 4^d cell on its own, so single cells scattered over the input compress in the sample to exactly
     * the bytes they take in the whole input. Its accuracy mode only looks at the binary exponent of the tolerance,
     * and its errors stay far below the tolerance.
     */
    TuningParams tuning_defaults(size_t /*ndims*/) const override {
        TuningParams p;
        p.independent_cell = 4;
        p.pow2_error_bounds = 1;
        p.rare_outlier_share = 0.01;
        return p;
    }
};

}  // namespace

const Backend &zfp_backend() {
    static const ZFPBackend b;
    return b;
}

}  // namespace qpet
