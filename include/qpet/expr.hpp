#ifndef QPET_EXPR_HPP
#define QPET_EXPR_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace qpet {

/**
 * Value and first two derivatives of a univariate function at one point.
 * dist estimates the distance (in x) from the point to the nearest singularity of the function, and brk the distance
 * to the nearest point where a condition of a piecewise function flips; +inf if there is none.
 */
struct Jet {
    double v = 0, d1 = 0, d2 = 0;
    double dist = std::numeric_limits<double>::infinity();
    double brk = std::numeric_limits<double>::infinity();
};

/**
 * A univariate QoI function f(x) given as text, e.g. "x^2", "log2(x)", "sin(10x)", "1/(1+exp(-x))".
 *
 * Grammar: + - * / ^ (or **), unary minus, parentheses, |...| for abs, implicit multiplication after a number ("10x"),
 * constants pi and e, and the functions sin cos tan asin acos atan sinh cosh tanh exp log ln log2 log10 sqrt cbrt abs
 * sigmoid, log(x, base), pow(a, b), min(a, b) and max(a, b), and the step functions floor ceil round sign.
 *
 * Piecewise functions: "condition ? a : b" or if(condition, a, b), with the comparisons < <= > >= == != and the words
 * and, or, not; e.g. "x < 0 ? 0 : sqrt(x)", "x < 1 ? x^2 : x < 5 ? x : 5".
 *
 * The text is compiled once to a postfix program. Derivatives come from second-order forward-mode
 * automatic differentiation, so no symbolic algebra library is needed.
 */
class Expr {
   public:
    /** @throws std::invalid_argument if the text does not parse */
    explicit Expr(const std::string &text);

    const std::string &text() const { return text_; }

    double eval(double x) const;

    /** out[i] = f(x[i]). x and out may be the same array. */
    void eval(const double *x, double *out, size_t n) const;

    Jet jet(double x) const;

    /** true if f(x) = a*x + b for constants a and b */
    bool is_linear() const { return linear_; }

   private:
    struct Op {
        uint8_t code;
        double c;
    };
    std::string text_;
    std::vector<Op> prog_;
    size_t depth_ = 0;
    bool linear_ = false;

    friend struct ExprCompiler;
};

/**
 * Estimate of the largest eb such that |x' - x| <= eb implies |f(x') - f(x)| <= tau, capped at eb_cap.
 * This is Theorem 5.1 / 5.2 of the QPET paper (second-order Taylor), followed by a check of f at x - eb and x + eb that
 * halves eb while either end is out of tolerance. It is an estimate, not a guarantee: the compressor still validates
 * every decompressed value.
 */
double pointwise_eb(const Expr &f, double x, double tau, double eb_cap);

}  // namespace qpet

#endif
