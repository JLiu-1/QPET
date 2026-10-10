#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "qpet/expr.hpp"

static int failures = 0;

#define CHECK(cond, ...) \
    do { \
        if (!(cond)) { \
            failures++; \
            std::printf("FAIL %s:%d: ", __FILE__, __LINE__); \
            std::printf(__VA_ARGS__); \
            std::printf("\n"); \
        } \
    } while (0)

static bool close(double a, double b, double rel = 1e-9) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    return std::fabs(a - b) <= rel * std::max({1.0, std::fabs(a), std::fabs(b)});
}

struct Case {
    std::string text;
    std::function<double(double)> ref;
    double lo, hi;  // test interval
};

int main() {
    const std::vector<Case> cases = {
        {"x", [](double x) { return x; }, -5, 5},
        {"x^2", [](double x) { return x * x; }, -5, 5},
        {"x**3", [](double x) { return x * x * x; }, -5, 5},
        {"x^2+x+1", [](double x) { return x * x + x + 1; }, -5, 5},
        {"-x^2", [](double x) { return -(x * x); }, -5, 5},
        {"2^x", [](double x) { return std::pow(2.0, x); }, -5, 5},
        {"e^x", [](double x) { return std::exp(x); }, -5, 5},
        {"exp(-x)", [](double x) { return std::exp(-x); }, -5, 5},
        {"log2(x)", [](double x) { return std::log2(x); }, 0.1, 50},
        {"ln(x)", [](double x) { return std::log(x); }, 0.1, 50},
        {"Log(x, 2)", [](double x) { return std::log(x) / std::log(2.0); }, 0.1, 50},
        {"log10(x)", [](double x) { return std::log10(x); }, 0.1, 50},
        {"1/x", [](double x) { return 1 / x; }, 0.2, 5},
        {"(x+3)^-1", [](double x) { return 1 / (x + 3); }, -2, 5},
        {"1/(x+1)", [](double x) { return 1 / (x + 1); }, 0, 5},
        {"sqrt(x)", [](double x) { return std::sqrt(x); }, 0.1, 50},
        {"x^0.5", [](double x) { return std::sqrt(x); }, 0.1, 50},
        {"sin(10x)", [](double x) { return std::sin(10 * x); }, -2, 2},
        {"sin(10*x)", [](double x) { return std::sin(10 * x); }, -2, 2},
        {"cos(x)*x", [](double x) { return std::cos(x) * x; }, -5, 5},
        {"tan(x)", [](double x) { return std::tan(x); }, -1, 1},
        {"tanh(x)", [](double x) { return std::tanh(x); }, -3, 3},
        {"(e^x-e^-x)/(e^x+e^-x)", [](double x) { return std::tanh(x); }, -3, 3},
        {"1/(1+exp(-x))", [](double x) { return 1 / (1 + std::exp(-x)); }, -5, 5},
        {"sigmoid(x)", [](double x) { return 1 / (1 + std::exp(-x)); }, -5, 5},
        {"|x|", [](double x) { return std::fabs(x); }, 0.5, 5},
        {"abs(x - 1)", [](double x) { return std::fabs(x - 1); }, 2, 5},
        {"pow(x, 2.5)", [](double x) { return std::pow(x, 2.5); }, 0.1, 5},
        {"x^x", [](double x) { return std::pow(x, x); }, 0.5, 3},
        {"atan(2x)", [](double x) { return std::atan(2 * x); }, -3, 3},
        {"3 - x/2 + pi", [](double x) { return 3 - x / 2 + M_PI; }, -5, 5},
        {"2(x+1)*x", [](double x) { return 2 * (x + 1) * x; }, -5, 5},
        {"cbrt(x)", [](double x) { return std::cbrt(x); }, 0.1, 5},
    };

    for (const Case &c : cases) {
        qpet::Expr f(c.text);
        const int steps = 97;
        std::vector<double> xs(steps), batch(steps);
        for (int i = 0; i < steps; i++) xs[i] = c.lo + (c.hi - c.lo) * (i + 0.37) / steps;
        f.eval(xs.data(), batch.data(), steps);
        for (int i = 0; i < steps; i++) {
            double x = xs[i], v = c.ref(x);
            CHECK(close(f.eval(x), v), "%s: eval(%g) = %.17g, want %.17g", c.text.c_str(), x, f.eval(x), v);
            CHECK(f.eval(x) == batch[i], "%s: batch and scalar eval differ at %g", c.text.c_str(), x);
            qpet::Jet j = f.jet(x);
            CHECK(close(j.v, v), "%s: jet value at %g", c.text.c_str(), x);
            double h = 1e-4 * std::max(1.0, std::fabs(x));
            double d1 = (c.ref(x + h) - c.ref(x - h)) / (2 * h);
            double d2 = (c.ref(x + h) - 2 * v + c.ref(x - h)) / (h * h);
            CHECK(close(j.d1, d1, 1e-5), "%s: f'(%g) = %.12g, finite difference %.12g", c.text.c_str(), x, j.d1, d1);
            CHECK(close(j.d2, d2, 2e-4), "%s: f''(%g) = %.12g, finite difference %.12g", c.text.c_str(), x, j.d2, d2);

            // the estimated bound holds on the whole interval for these functions
            for (double tau : {1e-1, 1e-3, 1e-6}) {
                double eb = qpet::pointwise_eb(f, x, tau, 1.0);
                CHECK(eb >= 0 && eb <= 1.0, "%s: eb(%g) = %g out of range", c.text.c_str(), x, eb);
                for (int k = -8; k <= 8; k++) {
                    double e = std::fabs(c.ref(x + eb * k / 8.0) - v);
                    // slack for the rounding of f at the two ends, where the error equals tau
                    CHECK(e <= tau * (1 + 1e-6) + 1e-12 * std::fabs(v), "%s: x=%g tau=%g eb=%g gives QoI error %g at %d/8", c.text.c_str(), x,
                          tau, eb, e, k);
                }
            }
        }
    }

    CHECK(qpet::Expr("x").is_linear(), "x is linear");
    CHECK(qpet::Expr("3 - x/2").is_linear(), "3 - x/2 is linear");
    CHECK(qpet::Expr("2(x+1)").is_linear(), "2(x+1) is linear");
    CHECK(!qpet::Expr("x^2").is_linear(), "x^2 is not linear");
    CHECK(!qpet::Expr("1/x").is_linear(), "1/x is not linear");

    // near a singularity the bound must not reach across it
    qpet::Expr lg("log2(x)");
    for (double x : {1e-6, 1e-3, 0.5}) {
        double eb = qpet::pointwise_eb(lg, x, 10.0, 1e9);
        CHECK(eb < x, "log2: eb(%g) = %g crosses 0", x, eb);
    }
    CHECK(qpet::pointwise_eb(lg, -1.0, 1.0, 1.0) == 0, "log2 of a negative value needs an exact value");
    CHECK(qpet::pointwise_eb(lg, 0.0, 1.0, 1.0) == 0, "log2 of 0 needs an exact value");

    for (const char *bad : {"", "x +", "foo(x)", "y", "(x", "2 ** ", "sin(x, 2)", "x x", "|x"}) {
        bool threw = false;
        try {
            qpet::Expr f(bad);
        } catch (const std::invalid_argument &) {
            threw = true;
        }
        CHECK(threw, "\"%s\" should not parse", bad);
    }

    std::printf(failures ? "%d FAILURES\n" : "all expression tests passed\n", failures);
    return failures ? 1 : 0;
}
