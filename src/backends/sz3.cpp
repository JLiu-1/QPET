#include <memory>
#include <stdexcept>

#include "SZ3/api/sz.hpp"
#include "qpet/backend.hpp"

namespace qpet {

namespace {

class SZ3Backend : public Backend {
   public:
    std::string name() const override { return "sz3"; }

    std::string version() const override { return SZ3_VER; }

    bool supports(DType, size_t ndims) const override { return ndims >= 1 && ndims <= 4; }

    /** option "algo": interp_lorenzo (default), interp, lorenzo_reg */
    std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims,
                                  const BackendParams &params) const override {
        SZ3::Config conf;
        conf.setDims(dims.begin(), dims.end());
        conf.errorBoundMode = SZ3::EB_ABS;
        conf.absErrorBound = params.abs_eb;
        conf.openmp = params.parallel && params.threads > 1;
        auto algo = params.options.find("algo");
        if (algo != params.options.end()) {
            if (algo->second == "interp_lorenzo") conf.cmprAlgo = SZ3::ALGO_INTERP_LORENZO;
            else if (algo->second == "interp") conf.cmprAlgo = SZ3::ALGO_INTERP;
            else if (algo->second == "lorenzo_reg") conf.cmprAlgo = SZ3::ALGO_LORENZO_REG;
            else throw std::invalid_argument("QPET: unknown sz3 algo " + algo->second);
        }
        size_t size = 0;
        std::unique_ptr<char[]> out(dtype == DType::F32
                                        ? SZ_compress<float>(conf, static_cast<const float *>(data), size)
                                        : SZ_compress<double>(conf, static_cast<const double *>(data), size));
        return std::vector<uint8_t>(out.get(), out.get() + size);
    }

    void decompress(const uint8_t *buf, size_t size, DType dtype, const std::vector<size_t> &, const BackendParams &,
                    void *out) const override {
        SZ3::Config conf;
        const char *in = reinterpret_cast<const char *>(buf);
        if (dtype == DType::F32) {
            float *p = static_cast<float *>(out);
            SZ_decompress<float>(conf, in, size, p);
        } else {
            double *p = static_cast<double *>(out);
            SZ_decompress<double>(conf, in, size, p);
        }
    }
};

}  // namespace

const Backend &sz3_backend() {
    static const SZ3Backend b;
    return b;
}

}  // namespace qpet
