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
        if (dims.size() == 2) {
            rc = C_API::sperr_comp_2d(data, is_float, dims[1], dims[0], MODE_PWE, params.abs_eb, 0, &dst, &len);
        } else {
            size_t chunk = 256;
            auto it = params.options.find("chunk");
            if (it != params.options.end()) chunk = std::stoul(it->second);
            rc = C_API::sperr_comp_3d(data, is_float, dims[2], dims[1], dims[0], chunk, chunk, chunk, MODE_PWE,
                                      params.abs_eb, static_cast<size_t>(params.threads), &dst, &len);
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
