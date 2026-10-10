#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <map>
#include <stdexcept>

#include <omp.h>
#include <zstd.h>

#include "qpet/qpet.hpp"

namespace qpet {

namespace {

const double INF = std::numeric_limits<double>::infinity();
const char MAGIC[4] = {'Q', 'P', 'E', 'T'};
const uint8_t FORMAT_VERSION = 2;
const size_t BATCH = 1024;

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

int resolve_threads(int t) { return t > 0 ? t : omp_get_max_threads(); }

// ---- byte streams -------------------------------------------------------------------------------------------------

struct Writer {
    std::vector<uint8_t> buf;

    template <class V>
    void put(V v) {
        const uint8_t *p = reinterpret_cast<const uint8_t *>(&v);
        buf.insert(buf.end(), p, p + sizeof(V));
    }
    void put_bytes(const void *p, size_t n) {
        const uint8_t *b = static_cast<const uint8_t *>(p);
        buf.insert(buf.end(), b, b + n);
    }
    void put_str(const std::string &s) {
        put<uint16_t>(static_cast<uint16_t>(s.size()));
        put_bytes(s.data(), s.size());
    }
};

struct Reader {
    const uint8_t *p, *end;

    const uint8_t *take(size_t n) {
        if (n > static_cast<size_t>(end - p)) throw std::runtime_error("QPET: compressed stream is truncated");
        const uint8_t *r = p;
        p += n;
        return r;
    }
    template <class V>
    V get() {
        V v;
        std::memcpy(&v, take(sizeof(V)), sizeof(V));
        return v;
    }
    std::string get_str() {
        size_t n = get<uint16_t>();
        return std::string(reinterpret_cast<const char *>(take(n)), n);
    }
};

void put_zstd(Writer &w, const void *data, size_t size) {
    std::vector<uint8_t> z(ZSTD_compressBound(size));
    size_t zs = ZSTD_compress(z.data(), z.size(), data, size, 3);
    if (ZSTD_isError(zs)) throw std::runtime_error("QPET: zstd compression failed");
    w.put<uint64_t>(size);
    w.put<uint64_t>(zs);
    w.put_bytes(z.data(), zs);
}

std::vector<uint8_t> get_zstd(Reader &r) {
    size_t size = r.get<uint64_t>(), zs = r.get<uint64_t>();
    const uint8_t *z = r.take(zs);
    if (ZSTD_getFrameContentSize(z, zs) != size) throw std::runtime_error("QPET: compressed stream is corrupted");
    std::vector<uint8_t> out(size);
    size_t got = ZSTD_decompress(out.data(), size, z, zs);
    if (ZSTD_isError(got) || got != size) throw std::runtime_error("QPET: compressed stream is corrupted");
    return out;
}

// ---- the values stored losslessly -----------------------------------------------------------------------------------

template <class T>
struct Patch {
    std::vector<uint64_t> idx;  // increasing
    std::vector<T> val;
};

template <class T>
void write_patch(Writer &w, const Patch<T> &patch) {
    w.put<uint64_t>(patch.idx.size());
    if (patch.idx.empty()) return;
    std::vector<uint8_t> gaps;  // LEB128 of the index differences
    gaps.reserve(patch.idx.size() * 2);
    uint64_t prev = 0;
    for (uint64_t i : patch.idx) {
        uint64_t d = i - prev;
        prev = i;
        while (d >= 0x80) {
            gaps.push_back(static_cast<uint8_t>(d) | 0x80);
            d >>= 7;
        }
        gaps.push_back(static_cast<uint8_t>(d));
    }
    put_zstd(w, gaps.data(), gaps.size());
    put_zstd(w, patch.val.data(), patch.val.size() * sizeof(T));
}

template <class T>
size_t apply_patch(Reader &r, T *out, size_t num) {
    size_t count = r.get<uint64_t>();
    if (count == 0) return 0;
    std::vector<uint8_t> gaps = get_zstd(r), vals = get_zstd(r);
    if (vals.size() != count * sizeof(T)) throw std::runtime_error("QPET: compressed stream is corrupted");
    size_t g = 0;
    uint64_t pos = 0;
    for (size_t k = 0; k < count; k++) {
        uint64_t d = 0;
        for (int shift = 0;; shift += 7) {
            if (g >= gaps.size() || shift > 63) throw std::runtime_error("QPET: compressed stream is corrupted");
            uint8_t b = gaps[g++];
            d |= static_cast<uint64_t>(b & 0x7F) << shift;
            if (!(b & 0x80)) break;
        }
        pos += d;
        if (pos >= num) throw std::runtime_error("QPET: compressed stream is corrupted");
        std::memcpy(out + pos, vals.data() + k * sizeof(T), sizeof(T));
    }
    return count;
}

// ---- validation -----------------------------------------------------------------------------------------------------

struct Errors {
    double qoi = 0, data = 0;
    void merge(const Errors &o) {
        qoi = std::max(qoi, o.qoi);
        data = std::max(data, o.data);
    }
};

/**
 * Checks dec[lo, hi) against orig. A value that breaks the QoI tolerance or the data bound is replaced by the original
 * and recorded in patch (with index base + i). Returns the largest errors among the values that were kept.
 */
template <class T>
Errors validate_range(const Expr &f, double tau, double eps, const T *orig, T *dec, size_t lo, size_t hi,
                      size_t base, Patch<T> &patch) {
    double xo[BATCH], xd[BATCH], qo[BATCH], qd[BATCH];
    Errors err;
    for (size_t s = lo; s < hi; s += BATCH) {
        size_t m = std::min(BATCH, hi - s);
        if (std::memcmp(orig + s, dec + s, m * sizeof(T)) == 0) continue;
        for (size_t i = 0; i < m; i++) {
            xo[i] = orig[s + i];
            xd[i] = dec[s + i];
        }
        f.eval(xo, qo, m);
        f.eval(xd, qd, m);
        for (size_t i = 0; i < m; i++) {
            if (xo[i] == xd[i] || (std::isnan(xo[i]) && std::isnan(xd[i]))) continue;
            double de = std::fabs(xo[i] - xd[i]), qe = std::fabs(qo[i] - qd[i]);
            // qe is NaN or inf when f is not finite at either value: such a value must be exact
            if (de <= eps && qe <= tau) {
                err.data = std::max(err.data, de);
                err.qoi = std::max(err.qoi, qe);
            } else {
                patch.idx.push_back(base + s + i);
                patch.val.push_back(orig[s + i]);
                dec[s + i] = orig[s + i];
            }
        }
    }
    return err;
}

template <class T>
Errors validate(const Expr &f, double tau, double eps, const T *orig, T *dec, size_t num, int threads,
                Patch<T> &patch) {
    threads = static_cast<int>(std::min<size_t>(threads, num / BATCH + 1));
    std::vector<Patch<T>> parts(threads);
    std::vector<Errors> errs(threads);
#pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (int t = 0; t < threads; t++) {
        size_t lo = num * t / threads, hi = num * (t + 1) / threads;
        errs[t] = validate_range(f, tau, eps, orig, dec, lo, hi, 0, parts[t]);
    }
    Errors err;
    for (int t = 0; t < threads; t++) {
        err.merge(errs[t]);
        patch.idx.insert(patch.idx.end(), parts[t].idx.begin(), parts[t].idx.end());
        patch.val.insert(patch.val.end(), parts[t].val.begin(), parts[t].val.end());
    }
    return err;
}

// ---- regional QoIs --------------------------------------------------------------------------------------------------

/**
 * The regions of a regional QoI: the array is cut into blocks of `edge` values per dimension, starting at the
 * origin. The blocks at the far ends are smaller when a dimension is not a multiple of the edge.
 */
struct Regions {
    std::vector<size_t> dims, counts, strides;
    size_t edge, total = 1;

    Regions(const std::vector<size_t> &dims_, size_t edge_) : dims(dims_), strides(dims_.size(), 1), edge(edge_) {
        for (size_t d : dims) {
            counts.push_back((d + edge - 1) / edge);
            total *= counts.back();
        }
        for (size_t d = dims.size() - 1; d > 0; d--) strides[d - 1] = strides[d] * dims[d];
    }

    /** flat indices of the values of region r, ascending */
    void members(size_t r, std::vector<size_t> &idx) const {
        const size_t N = dims.size();
        size_t origin[4], size[4], at[4] = {0, 0, 0, 0};
        for (size_t d = N; d-- > 0;) {
            origin[d] = (r % counts[d]) * edge;
            size[d] = std::min(edge, dims[d] - origin[d]);
            r /= counts[d];
        }
        idx.clear();
        for (;;) {
            size_t base = 0;
            for (size_t d = 0; d + 1 < N; d++) base += (origin[d] + at[d]) * strides[d];
            base += origin[N - 1];
            for (size_t i = 0; i < size[N - 1]; i++) idx.push_back(base + i);
            size_t d = N - 1;
            for (;;) {
                if (d == 0) return;
                d--;
                if (++at[d] < size[d]) break;
                at[d] = 0;
            }
        }
    }
};

/** per-thread scratch of the regional passes */
struct RegionScratch {
    std::vector<size_t> idx, order;
    std::vector<double> xo, xd, qo, qd, diff;
};

/**
 * The error of a regional QoI is |mean f(x') - mean f(x)| over the region. A value where f is not finite has to be
 * exact and takes no part in the means. This fills s.diff with f(x') - f(x) per value (0 for exact values) and
 * returns the sum; `bad` gets the values that must be restored whatever the mean is: f not finite at x or x', or
 * the data bound broken.
 */
template <class T>
double region_diffs(const Expr &f, double eps, const T *orig, const T *dec, RegionScratch &s,
                    std::vector<size_t> *bad, double *max_data_err) {
    const size_t n = s.idx.size();
    s.xo.resize(n), s.xd.resize(n), s.qo.resize(n), s.qd.resize(n), s.diff.resize(n);
    for (size_t i = 0; i < n; i++) {
        s.xo[i] = orig[s.idx[i]];
        s.xd[i] = dec[s.idx[i]];
    }
    f.eval(s.xo.data(), s.qo.data(), n);
    f.eval(s.xd.data(), s.qd.data(), n);
    double sum = 0;
    for (size_t i = 0; i < n; i++) {
        s.diff[i] = 0;
        if (s.xo[i] == s.xd[i] || (std::isnan(s.xo[i]) && std::isnan(s.xd[i]))) continue;
        double de = std::fabs(s.xo[i] - s.xd[i]), d = s.qd[i] - s.qo[i];
        if (de <= eps && std::isfinite(d)) {
            s.diff[i] = d;
            sum += d;
            if (max_data_err) *max_data_err = std::max(*max_data_err, de);
        } else if (bad) {
            bad->push_back(i);
        } else {  // verification: such a value is an unbounded error
            sum = INF;
            if (max_data_err) *max_data_err = std::max(*max_data_err, std::isnan(de) ? INF : de);
        }
    }
    return sum;
}

/**
 * Regional counterpart of validate(): every region whose mean QoI error is above tau gets values restored until it
 * is not. The values restored first are the ones that contribute most to the error, so that as few as possible are
 * stored losslessly.
 */
template <class T>
Errors validate_regions(const Expr &f, double tau, double eps, size_t edge, const T *orig, T *dec,
                        const std::vector<size_t> &dims, int threads, Patch<T> &patch) {
    const Regions regions(dims, edge);
    threads = static_cast<int>(std::min<size_t>(threads, regions.total));
    std::vector<std::vector<std::pair<size_t, T>>> restored(threads);
    std::vector<Errors> errs(threads);
#pragma omp parallel num_threads(threads)
    {
        RegionScratch s;
        std::vector<size_t> bad;
        Errors err;
        auto &mine = restored[omp_get_thread_num()];
        auto restore = [&](size_t i) {
            mine.emplace_back(s.idx[i], orig[s.idx[i]]);
            dec[s.idx[i]] = orig[s.idx[i]];
        };
#pragma omp for schedule(dynamic, 64)
        for (size_t r = 0; r < regions.total; r++) {
            regions.members(r, s.idx);
            const double n = static_cast<double>(s.idx.size());
            bad.clear();
            double data_err = 0;
            double sum = region_diffs(f, eps, orig, dec, s, &bad, &data_err);
            for (size_t i : bad) restore(i);
            // The loop ends on the sum as region_diffs() computes it, which is what verify() will compute too.
            while (std::fabs(sum) / n > tau) {
                s.order.clear();
                for (size_t i = 0; i < s.idx.size(); i++)
                    if (s.diff[i] != 0 && (s.diff[i] > 0) == (sum > 0)) s.order.push_back(i);
                std::sort(s.order.begin(), s.order.end(),
                          [&](size_t a, size_t b) { return std::fabs(s.diff[a]) > std::fabs(s.diff[b]); });
                double left = sum;
                for (size_t i : s.order) {
                    restore(i);
                    left -= s.diff[i];
                    if (std::fabs(left) / n <= tau) break;
                }
                data_err = 0;
                sum = region_diffs<T>(f, eps, orig, dec, s, nullptr, &data_err);
            }
            err.qoi = std::max(err.qoi, std::fabs(sum) / n);
            err.data = std::max(err.data, data_err);
        }
        errs[omp_get_thread_num()] = err;
    }
    Errors err;
    std::vector<std::pair<size_t, T>> all;
    for (int t = 0; t < threads; t++) {
        err.merge(errs[t]);
        all.insert(all.end(), restored[t].begin(), restored[t].end());
    }
    std::sort(all.begin(), all.end());
    for (const auto &v : all) {
        patch.idx.push_back(v.first);
        patch.val.push_back(v.second);
    }
    return err;
}

/** what has to hold for the decompressed data. region 1 is a point-wise QoI. */
struct Target {
    const Expr &f;
    double tau, eps;
    size_t region;
};

/** restores into dec, and records in patch, the values that have to be exact for the target to hold */
template <class T>
Errors enforce(const Target &t, const T *orig, T *dec, const std::vector<size_t> &dims, int threads,
               Patch<T> &patch) {
    if (t.region > 1) return validate_regions(t.f, t.tau, t.eps, t.region, orig, dec, dims, threads, patch);
    size_t num = 1;
    for (size_t d : dims) num *= d;
    return validate(t.f, t.tau, t.eps, orig, dec, num, threads, patch);
}

// ---- value ranges ---------------------------------------------------------------------------------------------------

struct Range {
    double lo = INF, hi = -INF;
    void add(double v) {
        if (std::isfinite(v)) {
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
    }
    void merge(const Range &o) {
        lo = std::min(lo, o.lo);
        hi = std::max(hi, o.hi);
    }
    double width() const { return hi >= lo ? hi - lo : 0; }
};

template <class T>
void ranges(const Expr &f, const T *data, size_t num, int threads, Range &data_range, Range &qoi_range) {
    std::vector<Range> dr(threads), qr(threads);
#pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (int t = 0; t < threads; t++) {
        double x[BATCH], q[BATCH];
        Range d, r;  // locals: writing dr[t] and qr[t] per value would share cache lines between threads
        size_t lo = num * t / threads, hi = num * (t + 1) / threads;
        for (size_t s = lo; s < hi; s += BATCH) {
            size_t m = std::min(BATCH, hi - s);
            for (size_t i = 0; i < m; i++) x[i] = data[s + i];
            f.eval(x, q, m);
            for (size_t i = 0; i < m; i++) {
                d.add(x[i]);
                r.add(q[i]);
            }
        }
        dr[t] = d;
        qr[t] = r;
    }
    for (int t = 0; t < threads; t++) {
        data_range.merge(dr[t]);
        qoi_range.merge(qr[t]);
    }
}

/** range of the regional QoI: the mean of f over each region, taken over the values where f is finite */
template <class T>
Range region_range(const Expr &f, const T *data, const std::vector<size_t> &dims, size_t edge, int threads) {
    const Regions regions(dims, edge);
    std::vector<Range> parts(threads);
#pragma omp parallel num_threads(threads)
    {
        std::vector<size_t> idx;
        std::vector<double> x, q;
        Range mine;
#pragma omp for schedule(static)
        for (size_t r = 0; r < regions.total; r++) {
            regions.members(r, idx);
            x.resize(idx.size()), q.resize(idx.size());
            for (size_t i = 0; i < idx.size(); i++) x[i] = data[idx[i]];
            f.eval(x.data(), q.data(), idx.size());
            double sum = 0;
            for (double v : q)
                if (std::isfinite(v)) sum += v;
            mine.add(sum / idx.size());
        }
        parts[omp_get_thread_num()] = mine;
    }
    Range range;
    for (const Range &p : parts) range.merge(p);
    return range;
}

// ---- tuning of the backend's error bound ----------------------------------------------------------------------------

/** point-wise error bounds of an evenly spaced sample of the data, sorted */
template <class T>
std::vector<double> sample_ebs(const Expr &f, const T *data, size_t num, double tau, double cap, int threads,
                               size_t samples) {
    size_t stride = std::max<size_t>(1, num / std::max<size_t>(samples, 1));
    size_t m = (num + stride - 1) / stride;
    std::vector<double> ebs(m);
#pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < m; i++) ebs[i] = pointwise_eb(f, data[i * stride], tau, cap);
    std::sort(ebs.begin(), ebs.end());
    return ebs;
}

double halton(size_t i, size_t base) {
    double r = 0, f = 1;
    for (i++; i > 0; i /= base) {
        f /= base;
        r += f * (i % base);
    }
    return r;
}

/** blocks of the input that stand in for it while candidate error bounds are tried */
template <class T>
struct Sample {
    std::vector<size_t> dims;  // of each block
    std::vector<std::vector<T>> blocks;
};

/** user settings over the backend's defaults over the general defaults: no field of the result is unset */
TuningParams resolve_tuning(const TuningParams &user, const TuningParams &backend, size_t N) {
    static const size_t EDGE[] = {1 << 18, 512, 64, 24};
    static const size_t MIN_EDGE[] = {1 << 18, 128, 64, 24};
    static const size_t MIN_BLOCKS[] = {5, 12, 5, 5};
    TuningParams p;
    auto pick = [](auto u, auto b, auto unset, auto fallback) { return u != unset ? u : b != unset ? b : fallback; };
    p.sample_rate = pick(user.sample_rate, backend.sample_rate, 0.0, 0.01);
    p.block_edge = pick(user.block_edge, backend.block_edge, size_t(0), EDGE[N - 1]);
    p.min_edge = pick(user.min_edge, backend.min_edge, size_t(0), std::min(p.block_edge, MIN_EDGE[N - 1]));
    p.min_blocks = pick(user.min_blocks, backend.min_blocks, size_t(0), MIN_BLOCKS[N - 1]);
    p.max_blocks = std::max(p.min_blocks, pick(user.max_blocks, backend.max_blocks, size_t(0), size_t(64)));
    p.independent_cell = pick(user.independent_cell, backend.independent_cell, size_t(0), size_t(1));
    p.pow2_error_bounds = pick(user.pow2_error_bounds, backend.pow2_error_bounds, -1, 0);
    p.refine_rounds = pick(user.refine_rounds, backend.refine_rounds, -1, 4);
    p.subtract_block_overhead = pick(user.subtract_block_overhead, backend.subtract_block_overhead, -1, 0);
    p.rare_outlier_share = pick(user.rare_outlier_share, backend.rare_outlier_share, -1.0, 0.1);
    p.quantiles = !user.quantiles.empty() ? user.quantiles : backend.quantiles;
    if (p.quantiles.empty()) p.quantiles = {0, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0};
    std::sort(p.quantiles.begin(), p.quantiles.end());
    p.point_samples = pick(user.point_samples, backend.point_samples, size_t(0), size_t(1) << 20);
    if (!(p.sample_rate > 0 && p.sample_rate <= 1)) throw std::invalid_argument("QPET: sample_rate must be in (0, 1]");
    for (double q : p.quantiles)
        if (!(q >= 0 && q <= 1)) throw std::invalid_argument("QPET: quantiles must be in [0, 1]");
    return p;
}

// copies the block of the given dimensions at `origin` (an offset into data) into out
template <class T>
void copy_block(const T *data, const std::vector<size_t> &strides, size_t origin, const std::vector<size_t> &bdims,
                T *out) {
    const size_t N = bdims.size(), row = bdims[N - 1];
    size_t rows = 1;
    for (size_t d = 0; d + 1 < N; d++) rows *= bdims[d];
    std::vector<size_t> at(N, 0);  // odometer over all dimensions but the last
    for (size_t r = 0; r < rows; r++) {
        size_t off = origin;
        for (size_t d = 0; d + 1 < N; d++) off += at[d] * strides[d];
        std::copy(data + off, data + off + row, out + r * row);
        for (size_t d = N - 1; d-- > 0;) {
            if (++at[d] < bdims[d]) break;
            at[d] = 0;
        }
    }
}

/**
 * The k-th block origin, as an offset into the data. Blocks sit on a tiling of each dimension (the last tile is
 * shifted back to fit) and the tiles are picked by a low-discrepancy sequence, so that every part of the input is
 * equally likely to be sampled. Uniformly random origins would favour the middle.
 */
size_t block_origin(size_t k, const std::vector<size_t> &dims, const std::vector<size_t> &bdims,
                    const std::vector<size_t> &strides, size_t align) {
    static const size_t PRIMES[] = {2, 3, 5, 7};
    size_t origin = 0;
    for (size_t d = 0; d < dims.size(); d++) {
        size_t tiles = (dims[d] + bdims[d] - 1) / bdims[d];
        size_t tile = std::min(tiles - 1, static_cast<size_t>(halton(k, PRIMES[d]) * tiles));
        size_t last = (dims[d] - bdims[d]) / align * align;  // origins stay multiples of align
        origin += std::min(tile * bdims[d], last) * strides[d];
    }
    return origin;
}

size_t lcm(size_t a, size_t b) {
    size_t x = a, y = b;
    while (y) {
        size_t t = x % y;
        x = y;
        y = t;
    }
    return a / x * b;
}

/**
 * align: edge of the regions of a regional QoI (1 for a point-wise one). Blocks then start on a region boundary and
 * hold whole regions, so that the regions of a block are regions of the input.
 */
template <class T>
Sample<T> sample_blocks(const T *data, const std::vector<size_t> &dims, const TuningParams &h, size_t align) {
    const double rate = h.sample_rate;
    const size_t N = dims.size();
    size_t num = 1, smallest = dims[0];
    for (size_t d : dims) {
        num *= d;
        smallest = std::min(smallest, d);
    }
    std::vector<size_t> strides(N, 1);
    for (size_t d = N - 1; d > 0; d--) strides[d - 1] = strides[d] * dims[d];
    Sample<T> s;

    const size_t cell_edge = h.independent_cell > 1 ? lcm(h.independent_cell, align) : 0;
    if (cell_edge && cell_edge <= 32 && smallest >= cell_edge) {
        // single cells, laid along the first dimension of one array
        const size_t c = cell_edge;
        std::vector<size_t> cell(N, c);
        size_t cell_num = 1;
        for (size_t d = 0; d < N; d++) cell_num *= c;
        size_t count = static_cast<size_t>(std::ceil(rate * num / cell_num));
        count = std::max<size_t>(count, 256);
        if (count * cell_num * 2 < num) {
            s.dims = cell;
            s.dims[0] = count * c;
            s.blocks.emplace_back(count * cell_num);
            for (size_t k = 0; k < count; k++)
                copy_block(data, strides, block_origin(k, dims, cell, strides, c), cell,
                           s.blocks[0].data() + k * cell_num);
            return s;
        }
    } else {
        // At least min_blocks blocks are taken, so on a small input blocks of the full edge would be a large share
        // of it. Shrink them until the sample is about `rate` of the input, down to min_edge.
        size_t fit = static_cast<size_t>(std::pow(rate * num / h.min_blocks, 1.0 / N));
        size_t edge = std::min(h.block_edge, std::max(fit, h.min_edge));
        size_t block_num = 1;
        for (size_t d : dims) {
            // the largest edge up to `edge` that splits the dimension into equal tiles
            size_t tiles = (d + edge - 1) / edge;
            size_t b = (d + tiles - 1) / tiles;
            if (align > 1 && d >= align) b = std::max(align, b / align * align);
            s.dims.push_back(b);
            block_num *= b;
        }
        size_t count = static_cast<size_t>(std::ceil(rate * num / block_num));
        count = std::min(std::max(count, h.min_blocks), h.max_blocks);
        if (count * block_num * 2 < num) {
            for (size_t k = 0; k < count; k++) {
                s.blocks.emplace_back(block_num);
                copy_block(data, strides, block_origin(k, dims, s.dims, strides, align), s.dims,
                           s.blocks[k].data());
            }
            return s;
        }
    }
    // small input: use all of it
    s.dims = dims;
    s.blocks.assign(1, std::vector<T>(data, data + num));
    return s;
}

template <class T>
class Tuner {
   public:
    /**
     * ebs: sorted point-wise error bounds of values spread evenly over the whole input.
     * point_tau: the point-wise QoI tolerance those bounds were computed for.
     */
    Tuner(const Backend &backend, DType dtype, const Sample<T> &sample, const Target &target, double point_tau,
          double cap, const std::vector<double> &ebs, const BackendParams &params, int threads,
          const TuningParams &hints)
        : backend_(backend), dtype_(dtype), sample_(sample), target_(target), ebs_(ebs), params_(params),
          threads_(threads), hints_(hints) {
        const Expr &f = target.f;
        const double tau = point_tau;
        params_.threads = 1;
        size_t num = sample.blocks[0].size();
        size_t stride = std::max<size_t>(1, sample.blocks.size() * num >> 20);
        size_t per_block = (num + stride - 1) / stride;
        block_ebs_.resize(sample.blocks.size() * per_block);
#pragma omp parallel for num_threads(threads) collapse(2)
        for (size_t b = 0; b < sample.blocks.size(); b++)
            for (size_t i = 0; i < per_block; i++)
                block_ebs_[b * per_block + i] = pointwise_eb(f, sample.blocks[b][i * stride], tau, cap);
        std::sort(block_ebs_.begin(), block_ebs_.end());
    }

    /**
     * Bytes per value the input is expected to take at this error bound: backend stream plus losslessly stored
     * values, both measured on the sample.
     *
     * The values that end up stored losslessly are the ones whose own error bound is below eb, and those can be
     * concentrated in a part of the input that the blocks over- or under-represent. So the bytes of the losslessly
     * stored values are scaled by how common such values are in the whole input relative to the blocks, and
     * estimated from the whole input alone if the blocks missed them.
     */
    double cost(double eb) {
        auto it = cache_.find(eb);
        if (it != cache_.end()) return it->second;
        BackendParams params = params_;
        params.abs_eb = eb;
        const size_t nb = sample_.blocks.size();
        size_t num = 1;
        for (size_t d : sample_.dims) num *= d;
        std::vector<size_t> base_bytes(nb, 0), patch_bytes(nb, 0), outliers(nb, 0);
        std::exception_ptr error;
#pragma omp parallel for num_threads(threads_) schedule(dynamic)
        for (size_t b = 0; b < nb; b++) {
            try {
                const T *orig = sample_.blocks[b].data();
                std::vector<uint8_t> cmp = backend_.compress(orig, dtype_, sample_.dims, params);
                std::vector<T> dec(num);
                backend_.decompress(cmp.data(), cmp.size(), dtype_, sample_.dims, params, dec.data());
                Patch<T> patch;
                enforce(target_, orig, dec.data(), sample_.dims, 1, patch);
                Writer w;
                write_patch(w, patch);
                base_bytes[b] = cmp.size();
                patch_bytes[b] = patch.idx.empty() ? 0 : w.buf.size();
                outliers[b] = patch.idx.size();
            } catch (...) {
#pragma omp critical
                error = std::current_exception();
            }
        }
        if (error) std::rethrow_exception(error);
        size_t overhead = 0;
        if (hints_.subtract_block_overhead && nb > 1) {
            std::vector<T> flat(num, sample_.blocks[0][0]);
            overhead = backend_.compress(flat.data(), dtype_, sample_.dims, params).size();
        }
        double base = 0, patch = 0;
        for (size_t b = 0; b < nb; b++) {
            base += base_bytes[b] > overhead ? base_bytes[b] - overhead : 0;
            patch += patch_bytes[b];
        }
        double in_blocks = below(block_ebs_, eb), in_input = below(ebs_, eb);
        double weight = 1;
        if (in_blocks > 0 && in_input <= 20 * in_blocks) {
            weight = std::max(0.05, in_input / in_blocks);
        } else if (in_input > 0) {
            // The blocks hold (almost) none of the values that cannot take this error bound, although the input
            // has them: rare values, clustered somewhere the blocks missed. What the blocks measured says nothing
            // about them, so a fixed share of them is assumed to end up stored losslessly.
            patch = in_input * hints_.rare_outlier_share * (sizeof(T) + 1) * nb * num;
        }
        double total = (base + patch * weight) / (nb * num);
        cache_[eb] = total;
        if (verbose)
        {
            size_t out = 0;
            for (size_t v : outliers) out += v;
            std::fprintf(stderr,
                         "QPET tuning: eb %-12g  backend %.4f + lossless %.4f x %.2f = %.4f bytes/value"
                         "  (stored losslessly: %.4f%% of the sample, sensitive: %.4f%% of the sample, %.4f%% of the input)\n",
                         eb, base / (nb * num), patch / (nb * num), weight, total, 100.0 * out / (nb * num),
                         100 * in_blocks, 100 * in_input);
        }
        return total;
    }

    int tried() const { return static_cast<int>(cache_.size()); }

    bool verbose = false;

   private:
    // fraction of the sorted values below v
    static double below(const std::vector<double> &sorted, double v) {
        return static_cast<double>(std::lower_bound(sorted.begin(), sorted.end(), v) - sorted.begin()) / sorted.size();
    }

    const Backend &backend_;
    DType dtype_;
    const Sample<T> &sample_;
    Target target_;
    const std::vector<double> &ebs_;
    std::vector<double> block_ebs_;
    BackendParams params_;
    int threads_;
    TuningParams hints_;
    std::map<double, double> cache_;
};

/**
 * Picks the backend's error bound. Candidates are quantiles of the point-wise error bounds: a small quantile leaves
 * few values to store losslessly, a large one lets the backend compress harder. Each candidate is tried on the sample
 * and the cheapest wins. The search then narrows around it, in log scale, because the bounds of some QoIs fall into
 * clusters far apart and the best choice can lie between two quantiles.
 */
template <class T>
double tune(Tuner<T> &tuner, const std::vector<double> &ebs, const TuningParams &hints) {
    std::vector<double> cand;
    for (double q : hints.quantiles) {
        double eb = ebs[std::min(ebs.size() - 1, static_cast<size_t>(q * ebs.size()))];
        if (eb <= 0) continue;
        if (hints.pow2_error_bounds) eb = std::ldexp(1.0, static_cast<int>(std::floor(std::log2(eb))));
        if (cand.empty() || eb > cand.back()) cand.push_back(eb);
    }
    if (cand.size() == 1 && !hints.pow2_error_bounds) return cand[0];
    size_t best = 0;
    for (size_t i = 1; i < cand.size(); i++)
        if (tuner.cost(cand[i]) < tuner.cost(cand[best])) best = i;
    double eb = cand[best];

    if (hints.pow2_error_bounds) {
        // Walk the exponent while it pays. The compressor may stay far inside its tolerance, so the best one can
        // lie above every point-wise bound: the walk is not stopped at the largest candidate.
        // The cost is not always unimodal in the exponent, so a direction is given up only after a few steps in a
        // row that do not pay.
        const double top = cand.back() * 4096, bottom = cand.front();
        for (double step : {2.0, 0.5}) {
            double probe = eb;
            for (int misses = 0; misses < 3;) {
                probe *= step;
                if (probe > top || probe < bottom) break;
                if (tuner.cost(probe) < tuner.cost(eb)) {
                    eb = probe;
                    misses = 0;
                } else {
                    misses++;
                }
            }
        }
        return eb;
    }

    double lo = best > 0 ? cand[best - 1] : eb, hi = best + 1 < cand.size() ? cand[best + 1] : eb;
    for (int round = 0; round < hints.refine_rounds && hi > lo * 1.05; round++) {
        double left = std::sqrt(lo * eb), right = std::sqrt(eb * hi);
        bool go_left = lo < eb && tuner.cost(left) < tuner.cost(eb);
        bool go_right = eb < hi && tuner.cost(right) < tuner.cost(go_left ? left : eb);
        if (go_right) {
            lo = eb;
            eb = right;
        } else if (go_left) {
            hi = eb;
            eb = left;
        } else {
            lo = left;
            hi = right;
        }
    }
    return eb;
}

// ---- container ------------------------------------------------------------------------------------------------------
//
// "QPET" | u8 format | u8 dtype | u8 ndims | u8 reserved | u64 dims[ndims]
// u64 region                                           (format 2 and later; format 1 is point-wise)
// f64 backend_eb | f64 qoi_tol_abs | f64 data_eb_abs
// str backend | str backend_version | str qoi          (str = u16 length + bytes)
// u64 base_bytes | backend stream
// u64 outliers | [zstd index gaps | zstd values]       (zstd = u64 raw size + u64 size + bytes)

struct Header {
    StreamInfo info;
    Reader body{nullptr, nullptr};  // positioned at base_bytes
};

Header read_header(const uint8_t *buf, size_t size) {
    Reader r{buf, buf + size};
    if (size < 8 || std::memcmp(r.take(4), MAGIC, 4) != 0) throw std::runtime_error("QPET: not a QPET stream");
    const uint8_t format = r.get<uint8_t>();
    if (format < 1 || format > FORMAT_VERSION) throw std::runtime_error("QPET: unsupported stream format version");
    Header h;
    uint8_t dtype = r.get<uint8_t>(), ndims = r.get<uint8_t>();
    r.get<uint8_t>();
    if (dtype > 1 || ndims < 1 || ndims > 4) throw std::runtime_error("QPET: compressed stream is corrupted");
    h.info.dtype = static_cast<DType>(dtype);
    h.info.num = 1;
    for (int d = 0; d < ndims; d++) {
        h.info.dims.push_back(r.get<uint64_t>());
        h.info.num *= h.info.dims.back();
    }
    h.info.region = format >= 2 ? r.get<uint64_t>() : 1;
    h.info.backend_eb = r.get<double>();
    h.info.qoi_tol_abs = r.get<double>();
    h.info.data_eb_abs = r.get<double>();
    h.info.backend = r.get_str();
    h.info.backend_version = r.get_str();
    h.info.qoi = r.get_str();
    h.body = r;
    return h;
}

template <class T>
std::vector<uint8_t> compress_impl(const T *data, DType dtype, const std::vector<size_t> &dims, const Config &conf,
                                   CompressInfo &info) {
    auto t_start = Clock::now();
    const Backend &backend = get_backend(conf.backend);
    if (dims.empty() || dims.size() > 4) throw std::invalid_argument("QPET: data must have 1 to 4 dimensions");
    size_t num = 1;
    for (size_t d : dims) num *= d;
    if (num == 0) throw std::invalid_argument("QPET: data is empty");
    if (!backend.supports(dtype, dims.size()))
        throw std::invalid_argument("QPET: backend " + conf.backend + " does not support this data type or " +
                                    std::to_string(dims.size()) + "-D data");
    if (!(conf.qoi_tol >= 0) || !(conf.data_eb >= 0))
        throw std::invalid_argument("QPET: tolerances must not be negative");

    const Expr f(conf.qoi);
    const int threads = resolve_threads(conf.threads);
    const size_t region = std::max<size_t>(conf.region, 1);

    // absolute tolerances
    auto t0 = Clock::now();
    Range data_range, qoi_range;
    ranges(f, data, num, threads, data_range, qoi_range);
    if (region > 1) qoi_range = region_range(f, data, dims, region, threads);
    double tau = conf.qoi_tol;
    if (conf.qoi_tol_mode == BoundMode::REL) tau *= qoi_range.width() > 0 ? qoi_range.width() : 1.0;
    double eps = INF;
    if (conf.data_eb > 0) {
        eps = conf.data_eb;
        if (conf.data_eb_mode == BoundMode::REL) eps *= data_range.width() > 0 ? data_range.width() : 1.0;
    }

    // The point-wise QoI tolerance that the estimates below start from. For the mean over a region of n values the
    // errors of the values partly cancel, so each may be larger than the tolerance of the mean: by the factor of
    // Theorem 5.4 of the paper (c = 2, confidence 0.999), and never less than the tolerance itself (Theorem 5.3).
    // This only seeds the candidates of the tuning; the guarantee comes from the check of every region.
    double point_tau = tau;
    if (region > 1) {
        double n = 1;
        for (size_t d : dims) n *= std::min(d, region);
        point_tau = tau * std::max(1.0, 2.0 * std::sqrt(n / (2 * std::log(2 / (1 - 0.999)))));
    }
    const Target target{f, tau, eps, region};

    // point-wise error bounds. A bound beyond the value range would not let the backend do anything more.
    double cap = std::min(eps, data_range.width());
    const TuningParams hints = resolve_tuning(conf.tuning, backend.tuning_defaults(dims.size()), dims.size());
    std::vector<double> ebs = sample_ebs(f, data, num, point_tau, cap, threads, hints.point_samples);
    info.t_estimate = since(t0);

    BackendParams params;
    params.threads = threads;
    params.parallel = conf.backend_parallel;
    params.options = conf.backend_options;

    t0 = Clock::now();
    double eb = 0;
    if (conf.backend_eb > 0) {
        eb = conf.backend_eb;
    } else if (ebs.back() <= 0) {
        // No value tolerates any error (or the data is constant): every lossy value will fail the check.
        eb = data_range.width() > 0 ? data_range.width() * 1e-9 : 1.0;
    } else if (conf.quantile >= 0) {
        size_t k = std::min(ebs.size() - 1, static_cast<size_t>(conf.quantile * ebs.size()));
        eb = ebs[k] > 0 ? ebs[k] : *std::upper_bound(ebs.begin(), ebs.end(), 0.0);
    } else {
        Sample<T> sample = sample_blocks(data, dims, hints, region);
        Tuner<T> tuner(backend, dtype, sample, target, point_tau, cap, ebs, params, threads, hints);
        tuner.verbose = conf.verbose;
        eb = tune(tuner, ebs, hints);
        info.candidates = tuner.tried();
    }
    info.t_tune = since(t0);
    params.abs_eb = eb;

    t0 = Clock::now();
    std::vector<uint8_t> base = backend.compress(data, dtype, dims, params);
    info.t_compress = since(t0);

    t0 = Clock::now();
    std::vector<T> dec(num);
    backend.decompress(base.data(), base.size(), dtype, dims, params, dec.data());
    info.t_decompress = since(t0);

    t0 = Clock::now();
    Patch<T> patch;
    Errors err = enforce(target, data, dec.data(), dims, threads, patch);
    info.t_validate = since(t0);

    Writer w;
    w.buf.reserve(base.size() + patch.idx.size() * (sizeof(T) + 2) + 256);
    w.put_bytes(MAGIC, 4);
    w.put<uint8_t>(FORMAT_VERSION);
    w.put<uint8_t>(static_cast<uint8_t>(dtype));
    w.put<uint8_t>(static_cast<uint8_t>(dims.size()));
    w.put<uint8_t>(0);
    for (size_t d : dims) w.put<uint64_t>(d);
    w.put<uint64_t>(region);
    w.put<double>(eb);
    w.put<double>(tau);
    w.put<double>(eps);
    w.put_str(backend.name());
    w.put_str(backend.version());
    w.put_str(f.text());
    w.put<uint64_t>(base.size());
    w.put_bytes(base.data(), base.size());
    size_t before_patch = w.buf.size();
    write_patch(w, patch);

    info.qoi_tol_abs = tau;
    info.data_eb_abs = eps;
    info.backend_eb = eb;
    info.quantile = static_cast<double>(std::lower_bound(ebs.begin(), ebs.end(), eb) - ebs.begin()) / ebs.size();
    info.num = num;
    info.outliers = patch.idx.size();
    info.base_bytes = base.size();
    info.patch_bytes = w.buf.size() - before_patch;
    info.total_bytes = w.buf.size();
    info.max_qoi_err = err.qoi;
    info.max_data_err = err.data;
    info.t_total = since(t_start);
    return std::move(w.buf);
}

}  // namespace

std::vector<uint8_t> compress(const void *data, DType dtype, const std::vector<size_t> &dims, const Config &conf,
                              CompressInfo *info) {
    CompressInfo local;
    CompressInfo &i = info ? *info : local;
    i = CompressInfo();
    if (dtype == DType::F32) return compress_impl(static_cast<const float *>(data), dtype, dims, conf, i);
    return compress_impl(static_cast<const double *>(data), dtype, dims, conf, i);
}

StreamInfo inspect(const uint8_t *buf, size_t size) {
    Header h = read_header(buf, size);
    Reader r = h.body;
    r.take(r.get<uint64_t>());
    h.info.outliers = r.get<uint64_t>();
    return h.info;
}

void decompress(const uint8_t *buf, size_t size, void *out, int threads) {
    Header h = read_header(buf, size);
    const StreamInfo &s = h.info;
    const Backend &backend = get_backend(s.backend);
    Reader r = h.body;
    size_t base_bytes = r.get<uint64_t>();
    const uint8_t *base = r.take(base_bytes);
    BackendParams params;
    params.abs_eb = s.backend_eb;
    params.threads = resolve_threads(threads);
    backend.decompress(base, base_bytes, s.dtype, s.dims, params, out);
    if (s.dtype == DType::F32) apply_patch(r, static_cast<float *>(out), s.num);
    else apply_patch(r, static_cast<double *>(out), s.num);
}

void TuningParams::set(const std::string &name, const std::string &value) {
    try {
        size_t used = 0;
        if (name == "quantiles") {
            quantiles.clear();
            for (size_t pos = 0; pos <= value.size();) {
                size_t end = std::min(value.find(':', pos), value.size());
                quantiles.push_back(std::stod(value.substr(pos, end - pos)));
                pos = end + 1;
            }
            return;
        }
        if (name == "sample_rate") sample_rate = std::stod(value, &used);
        else if (name == "block_edge") block_edge = std::stoull(value, &used);
        else if (name == "min_edge") min_edge = std::stoull(value, &used);
        else if (name == "min_blocks") min_blocks = std::stoull(value, &used);
        else if (name == "max_blocks") max_blocks = std::stoull(value, &used);
        else if (name == "independent_cell") independent_cell = std::stoull(value, &used);
        else if (name == "pow2_error_bounds") pow2_error_bounds = std::stoi(value, &used);
        else if (name == "refine_rounds") refine_rounds = std::stoi(value, &used);
        else if (name == "subtract_block_overhead") subtract_block_overhead = std::stoi(value, &used);
        else if (name == "point_samples") point_samples = std::stoull(value, &used);
        else if (name == "rare_outlier_share") rare_outlier_share = std::stod(value, &used);
        else throw std::invalid_argument("QPET: unknown tuning parameter \"" + name + "\"");
        if (used != value.size()) throw std::invalid_argument("");
    } catch (const std::invalid_argument &e) {
        if (std::string(e.what()).rfind("QPET", 0) == 0) throw;
        throw std::invalid_argument("QPET: bad value \"" + value + "\" for tuning parameter " + name);
    } catch (const std::out_of_range &) {
        throw std::invalid_argument("QPET: bad value \"" + value + "\" for tuning parameter " + name);
    }
}

void Config::set(const std::string &name, const std::string &value) {
    auto mode = [&]() {
        if (value == "abs" || value == "ABS") return BoundMode::ABS;
        if (value == "rel" || value == "REL") return BoundMode::REL;
        throw std::invalid_argument("QPET: " + name + " must be abs or rel, got \"" + value + "\"");
    };
    auto flag = [&]() {
        if (value == "0" || value == "false") return false;
        if (value == "1" || value == "true") return true;
        throw std::invalid_argument("QPET: " + name + " must be 0 or 1, got \"" + value + "\"");
    };
    if (name.rfind("tuning.", 0) == 0) return tuning.set(name.substr(7), value);
    if (name.rfind("backend.", 0) == 0) {
        backend_options[name.substr(8)] = value;
        return;
    }
    try {
        size_t used = value.size();
        if (name == "backend") backend = value;
        else if (name == "qoi") qoi = value;
        else if (name == "qoi_tol_mode") qoi_tol_mode = mode();
        else if (name == "data_eb_mode") data_eb_mode = mode();
        else if (name == "backend_parallel") backend_parallel = flag();
        else if (name == "verbose") verbose = flag();
        else if (name == "region") region = std::stoull(value, &used);
        else if (name == "qoi_tol") qoi_tol = std::stod(value, &used);
        else if (name == "data_eb") data_eb = std::stod(value, &used);
        else if (name == "threads") threads = std::stoi(value, &used);
        else if (name == "quantile") quantile = std::stod(value, &used);
        else if (name == "backend_eb") backend_eb = std::stod(value, &used);
        else throw std::invalid_argument("QPET: unknown option \"" + name + "\"");
        if (used != value.size()) throw std::invalid_argument("");
    } catch (const std::invalid_argument &e) {
        if (std::string(e.what()).rfind("QPET", 0) == 0) throw;
        throw std::invalid_argument("QPET: bad value \"" + value + "\" for option " + name);
    } catch (const std::out_of_range &) {
        throw std::invalid_argument("QPET: bad value \"" + value + "\" for option " + name);
    }
}

namespace {

template <class T>
VerifyResult verify_impl(const T *orig, const T *dec, size_t num, const Expr &f, int threads) {
    std::vector<Errors> errs(threads);
    std::vector<Range> dr(threads), qr(threads);
    std::vector<double> sq(threads, 0);
#pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (int t = 0; t < threads; t++) {
        double xo[BATCH], xd[BATCH], qo[BATCH], qd[BATCH];
        Errors e;
        Range d, r;
        double sum = 0;
        size_t lo = num * t / threads, hi = num * (t + 1) / threads;
        for (size_t s = lo; s < hi; s += BATCH) {
            size_t m = std::min(BATCH, hi - s);
            for (size_t i = 0; i < m; i++) {
                xo[i] = orig[s + i];
                xd[i] = dec[s + i];
            }
            f.eval(xo, qo, m);
            f.eval(xd, qd, m);
            for (size_t i = 0; i < m; i++) {
                d.add(xo[i]);
                r.add(qo[i]);
                if (xo[i] == xd[i] || (std::isnan(xo[i]) && std::isnan(xd[i]))) continue;
                double de = std::fabs(xo[i] - xd[i]), qe = std::fabs(qo[i] - qd[i]);
                // an error that is not a number counts as infinite
                e.data = std::max(e.data, std::isnan(de) ? INF : de);
                e.qoi = std::max(e.qoi, std::isnan(qe) ? INF : qe);
                if (std::isfinite(de)) sum += de * de;
            }
        }
        errs[t] = e;
        dr[t] = d;
        qr[t] = r;
        sq[t] = sum;
    }
    Errors err;
    Range data_range, qoi_range;
    double sum = 0;
    for (int t = 0; t < threads; t++) {
        err.merge(errs[t]);
        data_range.merge(dr[t]);
        qoi_range.merge(qr[t]);
        sum += sq[t];
    }
    VerifyResult v;
    v.max_qoi_err = err.qoi;
    v.max_data_err = err.data;
    v.qoi_range = qoi_range.width();
    v.data_range = data_range.width();
    v.psnr = sum > 0 ? 20 * std::log10(v.data_range) - 10 * std::log10(sum / num) : INF;
    return v;
}

/** the regional QoI error on top of what the point-wise pass found about the data */
template <class T>
void verify_regions(const T *orig, const T *dec, const std::vector<size_t> &dims, const Expr &f, size_t region,
                    int threads, VerifyResult &v) {
    const Regions regions(dims, region);
    std::vector<double> worst(threads, 0);
#pragma omp parallel num_threads(threads)
    {
        RegionScratch s;
        double mine = 0;
#pragma omp for schedule(static)
        for (size_t r = 0; r < regions.total; r++) {
            regions.members(r, s.idx);
            double sum = region_diffs<T>(f, INF, orig, dec, s, nullptr, nullptr);
            mine = std::max(mine, std::fabs(sum) / s.idx.size());
        }
        worst[omp_get_thread_num()] = mine;
    }
    v.max_qoi_err = *std::max_element(worst.begin(), worst.end());
    v.qoi_range = region_range(f, orig, dims, region, threads).width();
}

template <class T>
VerifyResult verify_typed(const T *orig, const T *dec, const std::vector<size_t> &dims, const Expr &f, size_t region,
                          int threads) {
    size_t num = 1;
    for (size_t d : dims) num *= d;
    VerifyResult v = verify_impl(orig, dec, num, f, threads);
    if (region > 1) verify_regions(orig, dec, dims, f, region, threads, v);
    return v;
}

}  // namespace

VerifyResult verify(const void *orig, const void *dec, DType dtype, const std::vector<size_t> &dims,
                    const std::string &qoi, size_t region, int threads) {
    if (dims.empty() || dims.size() > 4) throw std::invalid_argument("QPET: data must have 1 to 4 dimensions");
    const Expr f(qoi);
    threads = resolve_threads(threads);
    region = std::max<size_t>(region, 1);
    if (dtype == DType::F32)
        return verify_typed(static_cast<const float *>(orig), static_cast<const float *>(dec), dims, f, region, threads);
    return verify_typed(static_cast<const double *>(orig), static_cast<const double *>(dec), dims, f, region, threads);
}

VerifyResult verify(const void *orig, const void *dec, DType dtype, size_t num, const std::string &qoi, int threads) {
    return verify(orig, dec, dtype, std::vector<size_t>{num}, qoi, 1, threads);
}

}  // namespace qpet
