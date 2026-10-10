#include <memory>
#include <stdexcept>

#include "SZo/api/sz.hpp"
#include "qpet/backend.hpp"

// SZo and SZ3 cannot be included in one translation unit: both define versionInt() and versionStr() in the global
// namespace.

namespace qpet {

namespace {

class SZoBackend : public Backend {
   public:
    std::string name() const override { return "szo"; }

    std::string version() const override { return SZo_VER; }

    bool supports(DType, size_t ndims) const override { return ndims >= 1 && ndims <= 4; }

    /** option "algo": interp_lorenzo (default), interp, lorenzo_reg */
    std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims,
                                  const BackendParams &params) const override {
        SZo::Config conf;
        conf.setDims(dims.begin(), dims.end());
        conf.errorBoundMode = SZo::EB_ABS;
        conf.absErrorBound = params.abs_eb;
        conf.openmp = params.parallel && params.threads > 1;
        auto algo = params.options.find("algo");
        if (algo != params.options.end()) {
            if (algo->second == "interp_lorenzo") conf.cmprAlgo = SZo::ALGO_INTERP_LORENZO;
            else if (algo->second == "interp") conf.cmprAlgo = SZo::ALGO_INTERP;
            else if (algo->second == "lorenzo_reg") conf.cmprAlgo = SZo::ALGO_LORENZO_REG;
            else throw std::invalid_argument("QPET: unknown szo algo " + algo->second);
        }
        return dtype == DType::F32 ? run<float>(conf, data) : run<double>(conf, data);
    }

    void decompress(const uint8_t *buf, size_t size, DType dtype, const std::vector<size_t> &, const BackendParams &,
                    void *out) const override {
        SZo::Config conf;
        const char *in = reinterpret_cast<const char *>(buf);
        if (dtype == DType::F32) {
            float *p = static_cast<float *>(out);
            SZ_decompress<float>(conf, in, size, p);
        } else {
            double *p = static_cast<double *>(out);
            SZ_decompress<double>(conf, in, size, p);
        }
    }

   private:
    template <class T>
    static std::vector<uint8_t> run(const SZo::Config &conf, const void *data) {
        // SZo compresses in place and overwrites its input, so it gets a copy.
        const T *in = static_cast<const T *>(data);
        std::vector<T> copy(in, in + conf.num);
        size_t size = 0;
        std::unique_ptr<char[]> out(SZ_compress<T>(conf, copy.data(), size));
        return std::vector<uint8_t>(out.get(), out.get() + size);
    }
};

}  // namespace

const Backend &szo_backend() {
    static const SZoBackend b;
    return b;
}

}  // namespace qpet
