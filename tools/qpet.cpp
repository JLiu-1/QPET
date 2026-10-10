#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "qpet/qpet.hpp"

namespace {

void usage() {
    std::string backends;
    for (const std::string &b : qpet::available_backends()) backends += (backends.empty() ? "" : ", ") + b;
    std::printf(
        "QPET: QoI-preserving error-bounded lossy compression on top of an existing compressor\n\n"
        "Compress (and optionally decompress and verify):\n"
        "  qpet -b <backend> -f|-d -<N> <n1> ... <nN> -i <input> -q <QoI> -m <mode> <tol>\n"
        "       [-M <mode> <bound>] [-z <compressed>] [-o <decompressed>] [-a] [options]\n"
        "Decompress:\n"
        "  qpet -z <compressed> -o <decompressed>\n\n"
        "  -b <backend>        one of: %s\n"
        "  -f | -d             input is float32 | float64, raw binary\n"
        "  -1 nx | -2 nx ny | -3 nx ny nz | -4 nx ny nz nw\n"
        "                      dimensions, fastest-varying first (a 100x200x300 C array is -3 300 200 100)\n"
        "  -i <file>           input data\n"
        "  -z <file>           compressed file\n"
        "  -o <file>           decompressed data\n"
        "  -q <QoI>            point-wise QoI f(x), e.g. \"x^2\", \"log2(x)\", \"sin(10x)\", \"tanh(x)\"\n"
        "  -r <edge>           regional QoI: bound the error of the mean of f over blocks of <edge> values per\n"
        "                      dimension, instead of the error of f at every value\n"
        "  -m ABS|REL <tol>    QoI error threshold. REL is relative to the value range of the QoI.\n"
        "  -M ABS|REL <bound>  error bound on the data itself. REL is relative to the value range. Default: none.\n"
        "  -a                  decompress and report the errors\n"
        "  -t <n>              threads (default: all)\n"
        "  -P                  let SZ3/SZo run in their OpenMP mode (faster, lower compression ratio)\n"
        "  -Q <quantile>       skip the tuning: use this quantile of the point-wise error bounds\n"
        "  -E <bound>          skip the tuning: run the backend with this absolute error bound\n"
        "  -s <rate>           fraction of the data the tuning samples (default 0.01)\n"
        "  -S <name>=<value>   tuning parameter: sample_rate, block_edge, min_edge, min_blocks, max_blocks,\n"
        "                      independent_cell, pow2_error_bounds, refine_rounds, subtract_block_overhead,\n"
        "                      point_samples, rare_outlier_share, quantiles (e.g. -S quantiles=0:0.01:0.1:1). Unset ones follow the backend.\n"
        "  -O <key>=<value>    backend option, e.g. -O algo=interp (sz3, szo), -O chunk=128 (sperr)\n"
        "  -v                  print what the tuning tried\n"
        "  -l                  list the backends of this build\n",
        backends.c_str());
}

[[noreturn]] void die(const std::string &msg) {
    std::fprintf(stderr, "qpet: %s\n", msg.c_str());
    std::exit(1);
}

std::vector<uint8_t> read_file(const std::string &path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) die("cannot open " + path);
    std::vector<uint8_t> buf(static_cast<size_t>(in.tellg()));
    in.seekg(0);
    in.read(reinterpret_cast<char *>(buf.data()), buf.size());
    if (!in) die("cannot read " + path);
    return buf;
}

void write_file(const std::string &path, const void *data, size_t size) {
    std::ofstream out(path, std::ios::binary);
    out.write(static_cast<const char *>(data), size);
    if (!out) die("cannot write " + path);
}

qpet::BoundMode parse_mode(const std::string &s) {
    if (s == "ABS") return qpet::BoundMode::ABS;
    if (s == "REL") return qpet::BoundMode::REL;
    die("error mode must be ABS or REL, got " + s);
}

double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

}  // namespace

int main(int argc, char **argv) {
    qpet::Config conf;
    std::string input, compressed, output;
    std::vector<size_t> dims;  // as given: fastest first
    int dtype = -1;
    bool analyze = false, have_qoi_tol = false;

    if (argc == 1) {
        usage();
        return 0;
    }
    try {
        for (int i = 1; i < argc; i++) {
            std::string a = argv[i];
            auto next = [&]() -> std::string {
                if (i + 1 >= argc) die("option " + a + " needs a value");
                return argv[++i];
            };
            if (a == "-h" || a == "--help") {
                usage();
                return 0;
            } else if (a == "-l") {
                for (const std::string &b : qpet::available_backends())
                    std::printf("%s %s\n", b.c_str(), qpet::get_backend(b).version().c_str());
                return 0;
            } else if (a == "-b") conf.backend = next();
            else if (a == "-f") dtype = 0;
            else if (a == "-d") dtype = 1;
            else if (a == "-1" || a == "-2" || a == "-3" || a == "-4") {
                int n = a[1] - '0';
                for (int k = 0; k < n; k++) dims.push_back(std::stoull(next()));
            } else if (a == "-i") input = next();
            else if (a == "-z") compressed = next();
            else if (a == "-o") output = next();
            else if (a == "-q") conf.qoi = next();
            else if (a == "-m") {
                conf.qoi_tol_mode = parse_mode(next());
                conf.qoi_tol = std::stod(next());
                have_qoi_tol = true;
            } else if (a == "-M") {
                conf.data_eb_mode = parse_mode(next());
                conf.data_eb = std::stod(next());
            } else if (a == "-r") conf.region = std::stoull(next());
            else if (a == "-a") analyze = true;
            else if (a == "-t") conf.threads = std::stoi(next());
            else if (a == "-P") conf.backend_parallel = true;
            else if (a == "-Q") conf.quantile = std::stod(next());
            else if (a == "-v") conf.verbose = true;
            else if (a == "-E") conf.backend_eb = std::stod(next());
            else if (a == "-s") conf.tuning.sample_rate = std::stod(next());
            else if (a == "-S") {
                std::string kv = next();
                size_t eq = kv.find('=');
                if (eq == std::string::npos) die("-S takes name=value");
                conf.tuning.set(kv.substr(0, eq), kv.substr(eq + 1));
            }
            else if (a == "-O") {
                std::string kv = next();
                size_t eq = kv.find('=');
                if (eq == std::string::npos) die("-O takes key=value");
                conf.backend_options[kv.substr(0, eq)] = kv.substr(eq + 1);
            } else die("unknown option " + a);
        }

        if (input.empty()) {  // decompress only
            if (compressed.empty() || output.empty()) die("nothing to do: give -i to compress, or -z and -o");
            std::vector<uint8_t> cmp = read_file(compressed);
            qpet::StreamInfo s = qpet::inspect(cmp.data(), cmp.size());
            std::vector<uint8_t> dec(s.num * qpet::dtype_size(s.dtype));
            auto t = std::chrono::steady_clock::now();
            qpet::decompress(cmp.data(), cmp.size(), dec.data(), conf.threads);
            std::printf("backend = %s %s\nQoI = %s\ndecompression time = %.3f s\n", s.backend.c_str(),
                        s.backend_version.c_str(), s.qoi.c_str(), seconds_since(t));
            write_file(output, dec.data(), dec.size());
            return 0;
        }

        if (dtype < 0) die("give the data type with -f or -d");
        if (dims.empty()) die("give the dimensions with -1, -2, -3 or -4");
        if (!have_qoi_tol) die("give the QoI error threshold with -m");
        const qpet::DType dt = dtype == 0 ? qpet::DType::F32 : qpet::DType::F64;
        const size_t esize = qpet::dtype_size(dt);
        std::vector<size_t> cdims(dims.rbegin(), dims.rend());
        size_t num = 1;
        for (size_t d : dims) num *= d;

        std::vector<uint8_t> data = read_file(input);
        if (data.size() != num * esize)
            die(input + " has " + std::to_string(data.size()) + " bytes, the dimensions need " +
                std::to_string(num * esize));

        qpet::CompressInfo info;
        std::vector<uint8_t> cmp = qpet::compress(data.data(), dt, cdims, conf, &info);
        if (!compressed.empty()) write_file(compressed, cmp.data(), cmp.size());

        const double mb = num * esize / 1e6;
        std::printf("backend = %s %s\n", conf.backend.c_str(), qpet::get_backend(conf.backend).version().c_str());
        std::printf("QoI = %s\n", conf.qoi.c_str());
        if (conf.region > 1) std::printf("QoI region = %zu\n", conf.region);
        std::printf("QoI tolerance (abs) = %.6g\n", info.qoi_tol_abs);
        std::printf("data error bound (abs) = %.6g\n", info.data_eb_abs);
        std::printf("backend error bound = %.6g\n", info.backend_eb);
        std::printf("backend error bound quantile = %.4f\n", info.quantile);
        std::printf("tuning candidates = %d\n", info.candidates);
        std::printf("outliers = %zu (%.4f%%)\n", info.outliers, 100.0 * info.outliers / num);
        std::printf("backend bytes = %zu\n", info.base_bytes);
        std::printf("outlier bytes = %zu\n", info.patch_bytes);
        std::printf("compressed bytes = %zu\n", info.total_bytes);
        std::printf("compression ratio = %.4f\n", static_cast<double>(num * esize) / info.total_bytes);
        std::printf("bit rate = %.4f\n", 8.0 * info.total_bytes / num);
        std::printf("time estimate = %.3f s\n", info.t_estimate);
        std::printf("time tune = %.3f s\n", info.t_tune);
        std::printf("time backend compress = %.3f s\n", info.t_compress);
        std::printf("time backend decompress = %.3f s\n", info.t_decompress);
        std::printf("time validate = %.3f s\n", info.t_validate);
        std::printf("compression time = %.3f s\n", info.t_total);
        std::printf("compression speed = %.2f MB/s\n", mb / info.t_total);

        if (analyze || !output.empty()) {
            std::vector<uint8_t> dec(num * esize);
            auto t = std::chrono::steady_clock::now();
            qpet::decompress(cmp.data(), cmp.size(), dec.data(), conf.threads);
            double td = seconds_since(t);
            std::printf("decompression time = %.3f s\n", td);
            std::printf("decompression speed = %.2f MB/s\n", mb / td);
            if (!output.empty()) write_file(output, dec.data(), dec.size());
            if (analyze) {
                qpet::VerifyResult v =
                    qpet::verify(data.data(), dec.data(), dt, cdims, conf.qoi, conf.region, conf.threads);
                std::printf("max data error = %.6g (rel %.6g)\n", v.max_data_err,
                            v.data_range > 0 ? v.max_data_err / v.data_range : 0.0);
                std::printf("max QoI error = %.6g (rel %.6g)\n", v.max_qoi_err,
                            v.qoi_range > 0 ? v.max_qoi_err / v.qoi_range : 0.0);
                std::printf("PSNR = %.4f\n", v.psnr);
                bool ok = v.max_qoi_err <= info.qoi_tol_abs && v.max_data_err <= info.data_eb_abs;
                std::printf("bounds respected = %s\n", ok ? "yes" : "NO");
                if (!ok) return 2;
            }
        }
    } catch (const std::exception &e) {
        die(e.what());
    }
    return 0;
}
