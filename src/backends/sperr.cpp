#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "SPERR_C_API.h"
#include "SperrConfig.h"
#include "qpet/backend.hpp"

namespace qpet {

namespace {

const int MODE_PWE = 3;  // fixed point-wise error

/** point-wise error mode, 2-D and 3-D only */
class SPERRBackend : public Backend {
   public:
    std::string name() const override { return "sperr"; }

    std::string version() const override {
        return std::to_string(SPERR_VERSION_MAJOR) + "." + std::to_string(SPERR_VERSION_MINOR) + "." +
               std::to_string(SPERR_VERSION_PATCH);
    }

    bool supports(DType, size_t ndims) const override { return ndims == 2 || ndims == 3; }

    /** option "chunk": edge of the chunks a 3-D volume is split into (default 256) */
    std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims,
                                  const BackendParams &params) const override {
        const int is_float = dtype == DType::F32;
        void *dst = nullptr;
        size_t len = 0;
        int rc;
        // SPERR turns the wavelet coefficients into integers of at most 64 bits with a step of 1.5 tolerances, and
        // the errors it corrects one by one into integers with a step of one. A tolerance that is too small for
        // the magnitude of the data overflows the first, one that is too large overflows the second: SPERR then
        // fails, or crashes. The tolerance is therefore kept between the two limits. (QPET checks the result
        // against its own bounds anyway; the tolerance the stream was made with is not needed to decompress it.)
        // Data that leaves no room between the limits, or holds values that are not finite, is refused here.
        size_t num = 1;
        for (size_t d : dims) num *= d;
        double largest = 0;
        bool finite = true;
        for (size_t i = 0; i < num; i++) {
            double v = std::fabs(is_float ? static_cast<const float *>(data)[i] : static_cast<const double *>(data)[i]);
            if (!std::isfinite(v)) finite = false;
            else if (v > largest) largest = v;
        }
        if (!finite) throw std::runtime_error("QPET: sperr cannot compress data with NaN or infinite values");
        const double lowest = std::ldexp(largest, -48), highest = std::ldexp(1.0, 56);
        if (lowest > highest) {
            char text[32];
            std::snprintf(text, sizeof(text), "%g", largest);
            throw std::runtime_error(std::string("QPET: sperr cannot compress data with magnitudes above 2e31 (largest: ") +
                                     text + ")");
        }
        const double tolerance = std::min(std::max(params.abs_eb, lowest), highest);
        if (dims.size() == 2) {
            rc = C_API::sperr_comp_2d(data, is_float, dims[1], dims[0], MODE_PWE, tolerance, 0, &dst, &len);
        } else {
            size_t chunk = 256;
            auto it = params.options.find("chunk");
            if (it != params.options.end()) chunk = std::stoul(it->second);
            rc = C_API::sperr_comp_3d(data, is_float, dims[2], dims[1], dims[0], chunk, chunk, chunk, MODE_PWE,
                                      tolerance, static_cast<size_t>(params.threads), &dst, &len);
        }
        if (rc != 0) {
            std::free(dst);
            throw std::runtime_error("QPET: sperr compression failed with code " + std::to_string(rc));
        }
        const uint8_t *p = static_cast<const uint8_t *>(dst);
        std::vector<uint8_t> out(p, p + len);
        std::free(dst);
        return out;
    }

    void decompress(const uint8_t *buf, size_t size, DType dtype, const std::vector<size_t> &dims,
                    const BackendParams &params, void *out) const override {
        const int is_float = dtype == DType::F32;
        void *dst = nullptr;
        size_t num = 1;
        for (size_t d : dims) num *= d;
        int rc;
        if (dims.size() == 2) {
            rc = C_API::sperr_decomp_2d(buf, size, is_float, dims[1], dims[0], &dst);
        } else {
            size_t x = 0, y = 0, z = 0;
            rc = C_API::sperr_decomp_3d(buf, size, is_float, static_cast<size_t>(params.threads), &x, &y, &z, &dst);
            if (rc == 0 && x * y * z != num) rc = -1;
        }
        if (rc != 0) {
            std::free(dst);
            throw std::runtime_error("QPET: sperr decompression failed with code " + std::to_string(rc));
        }
        std::memcpy(out, dst, num * dtype_size(dtype));
        std::free(dst);
    }

    /** Most of SPERR's errors are a small fraction of the tolerance, so few of the sensitive values break a QoI. */
    TuningParams tuning_defaults(size_t /*ndims*/) const override {
        TuningParams p;
        p.rare_outlier_share = 0.02;
        return p;
    }
};

}  // namespace

const Backend &sperr_backend() {
    static const SPERRBackend b;
    return b;
}

}  // namespace qpet
