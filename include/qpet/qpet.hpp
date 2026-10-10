#ifndef QPET_HPP
#define QPET_HPP

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "qpet/backend.hpp"
#include "qpet/expr.hpp"

namespace qpet {

enum class BoundMode { ABS, REL };

struct Config {
    std::string backend = "sz3";

    /** point-wise QoI f(x), see Expr for the syntax */
    std::string qoi = "x";

    /**
     * 1: point-wise QoI, the error of f at every value is bounded.
     * r > 1: regional QoI. The array is cut into blocks of r values per dimension, starting at the origin (the
     * blocks at the far ends are smaller when a dimension is not a multiple of r), and the error of the mean of f over
     * every block is bounded.
     */
    size_t region = 1;

    /**
     * QoI error threshold. Point-wise: max |f(x) - f(x')|. Regional: max over the regions of
     * |mean f(x) - mean f(x')|. REL multiplies by the value range of the QoI over the input (of f, or of the
     * regional means of f).
     */
    double qoi_tol = 1e-3;
    BoundMode qoi_tol_mode = BoundMode::REL;

    /** Bound on the data itself: max |x - x'|. REL multiplies by (max x - min x). 0 means no bound on the data. */
    double data_eb = 0;
    BoundMode data_eb_mode = BoundMode::REL;

    /** threads for QPET's own passes and for backends whose stream does not depend on it. 0 = all available. */
    int threads = 0;

    /** let SZ3/SZo run their OpenMP mode, which is faster and lowers their compression ratio */
    bool backend_parallel = false;

    /** If >= 0, use this quantile of the point-wise error bounds as the backend's error bound and skip the tuning. */
    double quantile = -1;

    /** If > 0, run the backend with this absolute error bound and skip the tuning. */
    double backend_eb = 0;

    /** how the tuning samples and searches; the fields left unset follow the backend, see TuningParams */
    TuningParams tuning;

    std::map<std::string, std::string> backend_options;

    /** print what the tuning tried to stderr */
    bool verbose = false;

    /**
     * Sets a field by its name and a text value, e.g. set("qoi_tol_mode", "abs"), set("tuning.block_edge", "32"),
     * set("backend.algo", "interp") for a backend option.
     * @throws std::invalid_argument on an unknown name or a bad value
     */
    void set(const std::string &name, const std::string &value);
};

struct CompressInfo {
    double qoi_tol_abs = 0;
    double data_eb_abs = 0;  ///< +inf if there is no bound on the data
    double backend_eb = 0;   ///< absolute error bound the backend ran with
    double quantile = 0;     ///< where backend_eb falls among the point-wise error bounds
    int candidates = 0;      ///< error bounds the tuning tried
    size_t num = 0;
    size_t outliers = 0;     ///< values stored losslessly
    size_t base_bytes = 0;   ///< backend stream
    size_t patch_bytes = 0;  ///< losslessly stored values
    size_t total_bytes = 0;
    double max_qoi_err = 0;   ///< after correction
    double max_data_err = 0;  ///< after correction
    // seconds
    double t_estimate = 0, t_tune = 0, t_compress = 0, t_decompress = 0, t_validate = 0, t_total = 0;
};

struct StreamInfo {
    std::string backend, backend_version, qoi;
    DType dtype = DType::F32;
    std::vector<size_t> dims;
    size_t num = 0;
    size_t region = 1;  ///< 1: point-wise QoI
    double qoi_tol_abs = 0, data_eb_abs = 0, backend_eb = 0;
    size_t outliers = 0;
};

/**
 * QoI-preserving compression. The decompressed data x' satisfies |f(x) - f(x')| <= qoi_tol for every value (or, for
 * a regional QoI, |mean f(x) - mean f(x')| <= qoi_tol for every region) and, if a data bound is set,
 * |x - x'| <= data_eb for every value. A value where f is not finite is reproduced exactly. All of this holds
 * whatever the backend does, because the decompressed data is checked and the values that break a bound are stored
 * losslessly.
 * dims are in C order (slowest-varying first).
 */
std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims, const Config &conf,
                              CompressInfo *info = nullptr);

inline std::vector<uint8_t> compress(const float *data, const std::vector<size_t> &dims, const Config &conf,
                                     CompressInfo *info = nullptr) {
    return compress(data, DType::F32, dims, conf, info);
}

inline std::vector<uint8_t> compress(const double *data, const std::vector<size_t> &dims, const Config &conf,
                                     CompressInfo *info = nullptr) {
    return compress(data, DType::F64, dims, conf, info);
}

/** reads the header of a compressed stream */
StreamInfo inspect(const uint8_t *buf, size_t size);

/** out holds StreamInfo::num values of StreamInfo::dtype */
void decompress(const uint8_t *buf, size_t size, void *out, int threads = 0);

struct VerifyResult {
    double max_qoi_err = 0, max_data_err = 0;
    double qoi_range = 0, data_range = 0;
    double psnr = 0;  ///< of the data
};

/**
 * Errors between original and decompressed data, computed independently of compress().
 * With region > 1, max_qoi_err and qoi_range are those of the regional QoI (see Config::region).
 */
VerifyResult verify(const void *orig, const void *dec, DType dtype, const std::vector<size_t> &dims,
                    const std::string &qoi, size_t region = 1, int threads = 0);

/** point-wise QoI on a flat array */
VerifyResult verify(const void *orig, const void *dec, DType dtype, size_t num, const std::string &qoi,
                    int threads = 0);

}  // namespace qpet

#endif
