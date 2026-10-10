#ifndef QPET_BACKEND_HPP
#define QPET_BACKEND_HPP

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace qpet {

enum class DType : uint8_t { F32 = 0, F64 = 1 };

inline size_t dtype_size(DType t) { return t == DType::F32 ? 4 : 8; }

struct BackendParams {
    double abs_eb = 0;      ///< absolute error bound
    int threads = 1;        ///< threads the backend may use
    bool parallel = false;  ///< allow a parallel mode that changes the compressed stream (SZ3/SZo OpenMP)
    std::map<std::string, std::string> options;  ///< backend-specific, see each backend
};

/**
 * How the tuning of the backend's error bound samples the input and searches.
 *
 * The tuning compresses a sample instead of the whole input, and what a faithful sample looks like depends on how the
 * compressor works. Every field has an "unset" value (0, -1 or empty). A field is taken from, in this order: the
 * user's Config::tuning, the backend's Backend::tuning_defaults(), the defaults for the dimensionality.
 */
struct TuningParams {
    /** fraction of the input that is compressed for each candidate error bound (default 0.01) */
    double sample_rate = 0;

    /** edge of a sample block (default 1-D 262144, 2-D 512, 3-D 64, 4-D 24) */
    size_t block_edge = 0;

    /**
     * On a small input, min_blocks blocks of the full edge would be a large share of it. Blocks shrink towards this
     * edge until the sample is about sample_rate of the input (default: block_edge, except 128 in 2-D).
     */
    size_t min_edge = 0;

    /** bounds on the number of blocks (default 5 to 64, except 12 to 64 in 2-D) */
    size_t min_blocks = 0;
    size_t max_blocks = 0;

    /**
     * If > 1, the compressor codes cells of this edge independently of each other (zfp: 4). The sample is then a
     * large number of single cells laid side by side in one array: they compress to exactly the bytes they take in
     * the whole input, and they cover it far more evenly than a few blocks can. block_edge, min_edge and the block
     * counts are not used. 1 turns this off.
     */
    size_t independent_cell = 0;

    /** 1: only the binary exponent of the error bound changes the output (zfp), so try powers of two only. 0: off. */
    int pow2_error_bounds = -1;

    /** rounds of narrowing the search around the best candidate (default 4) */
    int refine_rounds = -1;

    /**
     * 1: subtract what an empty block costs from every sample block, for a compressor that writes a fixed amount per
     * call, which shows up once in the whole input but once per block in the sample. 0: off.
     */
    int subtract_block_overhead = -1;

    /**
     * Sample blocks can miss rare values that cannot take a candidate error bound, when those are clustered in a
     * small part of the input. The point samples still see them. This is the share of them assumed to end up stored
     * losslessly in that case; it depends on how far inside its error bound the compressor's errors stay
     * (default 0.1).
     */
    double rare_outlier_share = -1;

    /**
     * Quantiles of the point-wise error bounds that are the first candidates
     * (default 0, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1).
     */
    std::vector<double> quantiles;

    /** values whose point-wise error bound is computed to get those quantiles (default about 1 million) */
    size_t point_samples = 0;

    /**
     * Sets a field by its name, e.g. set("block_edge", "32") or set("quantiles", "0:0.01:0.1:1").
     * @throws std::invalid_argument on an unknown name or a bad value
     */
    void set(const std::string &name, const std::string &value);
};

/**
 * An error-bounded lossy compressor that QPET drives through its public API only.
 * Dimensions are in C order (slowest-varying first). Implementations must not modify the input and must be safe to
 * call from several threads at once.
 */
class Backend {
   public:
    virtual ~Backend() = default;

    virtual std::string name() const = 0;

    /** version of the compressor library this backend was built against */
    virtual std::string version() const = 0;

    virtual bool supports(DType dtype, size_t ndims) const = 0;

    virtual std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims,
                                          const BackendParams &params) const = 0;

    /** out holds prod(dims) values of dtype. params.abs_eb is the bound the stream was compressed with. */
    virtual void decompress(const uint8_t *buf, size_t size, DType dtype, const std::vector<size_t> &dims,
                            const BackendParams &params, void *out) const = 0;

    /** what suits this compressor; the fields it leaves unset take the general defaults */
    virtual TuningParams tuning_defaults(size_t /*ndims*/) const { return TuningParams(); }
};

/** @throws std::invalid_argument if this build has no backend of that name */
const Backend &get_backend(const std::string &name);

std::vector<std::string> available_backends();

}  // namespace qpet

#endif
