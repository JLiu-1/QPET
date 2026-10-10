#include "qpet/qpet_c.h"

#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>

#include "qpet/qpet.hpp"

namespace {

thread_local std::string last_error;

char *copy_text(const std::string &s) {
    char *p = static_cast<char *>(std::malloc(s.size() + 1));
    if (p) std::memcpy(p, s.c_str(), s.size() + 1);
    return p;
}

// runs f, turning an exception into an error code
template <class F>
int guarded(F f) {
    try {
        f();
        return 0;
    } catch (const std::exception &e) {
        last_error = e.what();
    } catch (...) {
        last_error = "QPET: unknown error";
    }
    return 1;
}

qpet::DType to_dtype(int dtype) {
    if (dtype != QPET_F32 && dtype != QPET_F64) throw std::invalid_argument("QPET: dtype must be QPET_F32 or QPET_F64");
    return dtype == QPET_F32 ? qpet::DType::F32 : qpet::DType::F64;
}

std::vector<size_t> to_dims(size_t ndims, const size_t *dims) {
    if (!dims || ndims < 1 || ndims > 4) throw std::invalid_argument("QPET: data must have 1 to 4 dimensions");
    return std::vector<size_t>(dims, dims + ndims);
}

// full precision, so that a tolerance read back compares exactly
std::string number(double v) {
    std::ostringstream o;
    o.precision(17);
    o << v;
    return o.str();
}

}  // namespace

extern "C" {

int qpet_compress(const void *data, int dtype, size_t ndims, const size_t *dims, const char *options,
                  unsigned char **out, size_t *out_size, char **info) {
    return guarded([&] {
        if (!data || !out || !out_size) throw std::invalid_argument("QPET: null argument");
        qpet::Config conf;
        std::istringstream lines(options ? options : "");
        for (std::string line; std::getline(lines, line);) {
            if (line.empty()) continue;
            size_t eq = line.find('=');
            if (eq == std::string::npos) throw std::invalid_argument("QPET: option \"" + line + "\" is not name=value");
            conf.set(line.substr(0, eq), line.substr(eq + 1));
        }
        qpet::CompressInfo i;
        std::vector<uint8_t> cmp = qpet::compress(data, to_dtype(dtype), to_dims(ndims, dims), conf, &i);
        unsigned char *buf = static_cast<unsigned char *>(std::malloc(cmp.size()));
        if (!buf) throw std::bad_alloc();
        std::memcpy(buf, cmp.data(), cmp.size());
        *out = buf;
        *out_size = cmp.size();
        if (info) {
            std::ostringstream o;
            o << "qoi_tol_abs=" << number(i.qoi_tol_abs) << "\ndata_eb_abs=" << number(i.data_eb_abs)
              << "\nbackend_eb=" << number(i.backend_eb) << "\nquantile=" << number(i.quantile)
              << "\ncandidates=" << i.candidates << "\nnum=" << i.num << "\noutliers=" << i.outliers
              << "\nbase_bytes=" << i.base_bytes << "\npatch_bytes=" << i.patch_bytes
              << "\ntotal_bytes=" << i.total_bytes << "\nmax_qoi_err=" << number(i.max_qoi_err)
              << "\nmax_data_err=" << number(i.max_data_err) << "\nt_estimate=" << number(i.t_estimate)
              << "\nt_tune=" << number(i.t_tune) << "\nt_compress=" << number(i.t_compress)
              << "\nt_decompress=" << number(i.t_decompress) << "\nt_validate=" << number(i.t_validate)
              << "\nt_total=" << number(i.t_total) << "\n";
            *info = copy_text(o.str());
        }
    });
}

int qpet_inspect(const unsigned char *buf, size_t size, char **info) {
    return guarded([&] {
        if (!buf || !info) throw std::invalid_argument("QPET: null argument");
        qpet::StreamInfo s = qpet::inspect(buf, size);
        std::ostringstream o;
        o << "backend=" << s.backend << "\nbackend_version=" << s.backend_version << "\nqoi=" << s.qoi
          << "\nregion=" << s.region << "\ndtype=" << (s.dtype == qpet::DType::F32 ? QPET_F32 : QPET_F64) << "\ndims=";
        for (size_t d = 0; d < s.dims.size(); d++) o << (d ? "," : "") << s.dims[d];
        o << "\nnum=" << s.num << "\nqoi_tol_abs=" << number(s.qoi_tol_abs) << "\ndata_eb_abs="
          << number(s.data_eb_abs) << "\nbackend_eb=" << number(s.backend_eb) << "\noutliers=" << s.outliers << "\n";
        *info = copy_text(o.str());
    });
}

int qpet_decompress(const unsigned char *buf, size_t size, void *out, size_t out_bytes, int threads) {
    return guarded([&] {
        if (!buf || !out) throw std::invalid_argument("QPET: null argument");
        qpet::StreamInfo s = qpet::inspect(buf, size);
        if (out_bytes != s.num * qpet::dtype_size(s.dtype))
            throw std::invalid_argument("QPET: the output buffer does not have the size of the decompressed data");
        qpet::decompress(buf, size, out, threads);
    });
}

int qpet_verify(const void *orig, const void *dec, int dtype, size_t ndims, const size_t *dims, const char *qoi,
                size_t region, int threads, char **result) {
    return guarded([&] {
        if (!orig || !dec || !qoi || !result) throw std::invalid_argument("QPET: null argument");
        qpet::VerifyResult v = qpet::verify(orig, dec, to_dtype(dtype), to_dims(ndims, dims), qoi, region, threads);
        std::ostringstream o;
        o << "max_qoi_err=" << number(v.max_qoi_err) << "\nmax_data_err=" << number(v.max_data_err)
          << "\nqoi_range=" << number(v.qoi_range) << "\ndata_range=" << number(v.data_range)
          << "\npsnr=" << number(v.psnr) << "\n";
        *result = copy_text(o.str());
    });
}

int qpet_backends(char **list) {
    return guarded([&] {
        if (!list) throw std::invalid_argument("QPET: null argument");
        std::string s;
        for (const std::string &b : qpet::available_backends()) s += b + "=" + qpet::get_backend(b).version() + "\n";
        *list = copy_text(s);
    });
}

const char *qpet_last_error(void) { return last_error.c_str(); }

void qpet_free(void *p) { std::free(p); }

}  // extern "C"
