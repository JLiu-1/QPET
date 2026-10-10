#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "qpet/qpet.hpp"

static int failures = 0;

// smooth waves plus a little noise, all values in about [1, 3]
template <class T>
std::vector<T> make_data(const std::vector<size_t> &dims, bool hostile) {
    size_t num = 1;
    for (size_t d : dims) num *= d;
    std::vector<T> data(num);
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> noise(-0.005, 0.005);
    std::vector<size_t> at(dims.size(), 0);
    for (size_t i = 0; i < num; i++) {
        double v = 2;
        for (size_t d = 0; d < dims.size(); d++) v += std::sin(0.07 * (d + 1) * at[d] + d) / dims.size();
        data[i] = static_cast<T>(v + noise(rng));
        for (size_t d = dims.size(); d-- > 0;) {
            if (++at[d] < dims[d]) break;
            at[d] = 0;
        }
    }
    if (hostile) {  // values where log2 is not finite, and a flat region
        for (size_t i = 0; i < num / 50; i++) data[i * 37 % num] = 0;
        for (size_t i = 0; i < num / 80; i++) data[i * 53 % num] = static_cast<T>(-1.5);
        for (size_t i = num / 2; i < num / 2 + num / 20; i++) data[i] = static_cast<T>(2.5);
    }
    return data;
}

// Largest error of the mean of f over the blocks of `region` values per dimension, written without any of the
// library's region code: a value where f is not finite must be exact and is left out of the mean.
template <class T>
double brute_region_error(const std::vector<T> &a, const std::vector<T> &b, const std::vector<size_t> &dims,
                          size_t region, const std::string &qoi) {
    qpet::Expr f(qoi);
    const size_t N = dims.size();
    std::vector<size_t> counts(N);
    size_t blocks = 1;
    for (size_t d = 0; d < N; d++) blocks *= counts[d] = (dims[d] + region - 1) / region;
    std::vector<double> sum(blocks, 0), cnt(blocks, 0);
    std::vector<size_t> at(N, 0);
    double worst = 0;
    for (size_t i = 0; i < a.size(); i++) {
        size_t blk = 0;
        for (size_t d = 0; d < N; d++) blk = blk * counts[d] + at[d] / region;
        double qa = f.eval(a[i]), qb = f.eval(b[i]);
        cnt[blk] += 1;
        if (a[i] != b[i]) {
            if (std::isfinite(qa) && std::isfinite(qb)) sum[blk] += qb - qa;
            else worst = INFINITY;
        }
        for (size_t d = N; d-- > 0;) {
            if (++at[d] < dims[d]) break;
            at[d] = 0;
        }
    }
    for (size_t k = 0; k < blocks; k++) worst = std::max(worst, std::fabs(sum[k]) / cnt[k]);
    return worst;
}

template <class T>
void run(const std::string &backend, const std::vector<size_t> &dims, const std::string &qoi, double tol, double eb,
         bool hostile, size_t region = 1) {
    const qpet::DType dtype = sizeof(T) == 4 ? qpet::DType::F32 : qpet::DType::F64;
    if (!qpet::get_backend(backend).supports(dtype, dims.size())) return;
    std::vector<T> data = make_data<T>(dims, hostile);
    qpet::Config conf;
    conf.backend = backend;
    conf.qoi = qoi;
    conf.qoi_tol = tol;
    conf.data_eb = eb;
    conf.region = region;
    conf.threads = 4;
    std::string what = backend + " " + (sizeof(T) == 4 ? "f32" : "f64") + " " + std::to_string(dims.size()) + "D " +
                       qoi + " tol=" + std::to_string(tol) + (region > 1 ? " region=" + std::to_string(region) : "") +
                       (hostile ? " hostile" : "");
    try {
        qpet::CompressInfo info;
        std::vector<uint8_t> cmp = qpet::compress(data.data(), dims, conf, &info);
        qpet::StreamInfo s = qpet::inspect(cmp.data(), cmp.size());
        std::vector<T> dec(data.size());
        qpet::decompress(cmp.data(), cmp.size(), dec.data());
        qpet::VerifyResult v = qpet::verify(data.data(), dec.data(), dtype, dims, qoi, region);
        bool ok = s.backend == backend && s.dims == dims && s.num == data.size() && s.outliers == info.outliers &&
                  s.region == region && v.max_qoi_err <= info.qoi_tol_abs && v.max_data_err <= info.data_eb_abs &&
                  cmp.size() == info.total_bytes;
        if (region > 1) {
            double brute = brute_region_error(data, dec, dims, region, qoi);
            // the sums are taken in a different order than in the library, hence the slack of a few ulps
            ok = ok && brute <= info.qoi_tol_abs * (1 + 1e-9) && std::fabs(brute - v.max_qoi_err) <= 1e-9 * brute + 1e-300;
        }
        std::printf("%s %-58s CR %7.2f  outliers %6.3f%%  QoI err/tol %.3f\n", ok ? "ok  " : "FAIL", what.c_str(),
                    static_cast<double>(data.size() * sizeof(T)) / cmp.size(), 100.0 * info.outliers / data.size(),
                    info.qoi_tol_abs > 0 ? v.max_qoi_err / info.qoi_tol_abs : 0.0);
        if (!ok) failures++;
    } catch (const std::exception &e) {
        std::printf("FAIL %s: %s\n", what.c_str(), e.what());
        failures++;
    }
}

int main() {
    const std::vector<std::string> qois = {"x",       "x^2",     "x^3",       "log2(x)", "sin(10x)",
                                           "tanh(x)", "1/x",     "sqrt(x)",   "exp(x)",  "x^2+log2(x)"};
    const std::vector<std::vector<size_t>> shapes = {{5000}, {120, 150}, {40, 56, 72}, {6, 20, 24, 30}, {300, 300, 300}};
    for (const std::string &backend : qpet::available_backends()) {
        for (const auto &dims : shapes) {
            bool big = dims[0] == 300;
            for (const std::string &qoi : qois) {
                if (big && qoi != "x^2" && qoi != "log2(x)") continue;
                run<float>(backend, dims, qoi, 1e-3, 1e-2, false);
                if (!big) run<double>(backend, dims, qoi, 1e-4, 0, false);
            }
            run<float>(backend, dims, "log2(x)", 1e-3, 1e-2, true);
            if (!big) run<double>(backend, dims, "sqrt(x)", 1e-3, 0, true);
        }
    }
    // regional QoIs: regions that divide the array, regions that do not, a region larger than a dimension
    for (const std::string &backend : qpet::available_backends()) {
        for (const auto &dims : std::vector<std::vector<size_t>>{{5000}, {120, 150}, {40, 56, 72}, {6, 20, 24, 30}}) {
            for (size_t region : {2, 4, 7, 16}) {
                for (const std::string &qoi : {"x", "x^2", "log2(x)", "tanh(x)"}) {
                    run<float>(backend, dims, qoi, 1e-3, 1e-2, false, region);
                    run<double>(backend, dims, qoi, 1e-4, 0, false, region);
                }
                run<float>(backend, dims, "log2(x)", 1e-3, 1e-2, true, region);
                run<double>(backend, dims, "sqrt(x)", 1e-3, 0, true, region);
            }
        }
        run<float>(backend, {300, 300, 300}, "x^2", 1e-3, 1e-2, false, 4);
        run<float>(backend, {40, 56, 72}, "x^2", 0, 0, false, 4);
    }

    // a tolerance of zero means lossless
    for (const std::string &backend : qpet::available_backends()) run<float>(backend, {40, 56, 72}, "x^2", 0, 0, false);

    std::printf(failures ? "%d FAILURES\n" : "all round trips passed\n", failures);
    return failures ? 1 : 0;
}
