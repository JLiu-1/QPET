#include "qpet/expr.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace qpet {

namespace {

const double INF = std::numeric_limits<double>::infinity();
const size_t MAX_DEPTH = 64;
const size_t CHUNK = 256;

enum OpCode : uint8_t {
    OP_CONST,
    OP_VAR,
    OP_ADD,
    OP_SUB,
    OP_MUL,
    OP_DIV,
    OP_POW,
    OP_MIN,
    OP_MAX,
    // conditions: the result is 1 (true) or 0
    OP_LT,
    OP_LE,
    OP_GT,
    OP_GE,
    OP_EQ,
    OP_NE,
    OP_AND,
    OP_OR,
    OP_SELECT,  // condition, value if true, value if false
    // one operand and a constant
    OP_ADDC,
    OP_MULC,
    OP_DIVC,   // u / c
    OP_RDIVC,  // c / u
    OP_POWC,   // u ^ c
    // one operand
    OP_NEG,
    OP_SIN,
    OP_COS,
    OP_TAN,
    OP_ASIN,
    OP_ACOS,
    OP_ATAN,
    OP_SINH,
    OP_COSH,
    OP_TANH,
    OP_EXP,
    OP_LOG,
    OP_LOG2,
    OP_LOG10,
    OP_SQRT,
    OP_CBRT,
    OP_ABS,
    OP_SIGMOID,
    OP_NOT,
    // steps: constant between breakpoints
    OP_FLOOR,
    OP_CEIL,
    OP_ROUND,
    OP_SIGN,
};

struct Func {
    const char *name;
    OpCode op;
};
const Func FUNCS[] = {
    {"sin", OP_SIN},     {"cos", OP_COS},   {"tan", OP_TAN},   {"asin", OP_ASIN}, {"acos", OP_ACOS},
    {"atan", OP_ATAN},   {"sinh", OP_SINH}, {"cosh", OP_COSH}, {"tanh", OP_TANH}, {"exp", OP_EXP},
    {"log", OP_LOG},     {"ln", OP_LOG},    {"log2", OP_LOG2}, {"log10", OP_LOG10}, {"sqrt", OP_SQRT},
    {"cbrt", OP_CBRT},   {"abs", OP_ABS},   {"sigmoid", OP_SIGMOID},
    {"floor", OP_FLOOR}, {"ceil", OP_CEIL}, {"round", OP_ROUND}, {"sign", OP_SIGN},
};

inline double sigmoid(double u) { return 1.0 / (1.0 + std::exp(-u)); }
inline double sign(double u) { return u > 0 ? 1.0 : u < 0 ? -1.0 : u; }  // 0 at 0, NaN stays NaN

inline double unary_value(OpCode op, double u) {
    switch (op) {
        case OP_NEG: return -u;
        case OP_SIN: return std::sin(u);
        case OP_COS: return std::cos(u);
        case OP_TAN: return std::tan(u);
        case OP_ASIN: return std::asin(u);
        case OP_ACOS: return std::acos(u);
        case OP_ATAN: return std::atan(u);
        case OP_SINH: return std::sinh(u);
        case OP_COSH: return std::cosh(u);
        case OP_TANH: return std::tanh(u);
        case OP_EXP: return std::exp(u);
        case OP_LOG: return std::log(u);
        case OP_LOG2: return std::log2(u);
        case OP_LOG10: return std::log10(u);
        case OP_SQRT: return std::sqrt(u);
        case OP_CBRT: return std::cbrt(u);
        case OP_ABS: return std::fabs(u);
        case OP_SIGMOID: return sigmoid(u);
        case OP_NOT: return u == 0;
        case OP_FLOOR: return std::floor(u);
        case OP_CEIL: return std::ceil(u);
        case OP_ROUND: return std::round(u);
        case OP_SIGN: return sign(u);
        default: return std::nan("");
    }
}

inline double binary_value(OpCode op, double a, double b) {
    switch (op) {
        case OP_ADD: return a + b;
        case OP_SUB: return a - b;
        case OP_MUL: return a * b;
        case OP_DIV: return a / b;
        case OP_POW: return std::pow(a, b);
        case OP_MIN: return a < b ? a : b;
        case OP_MAX: return a > b ? a : b;
        case OP_LT: return a < b;
        case OP_LE: return a <= b;
        case OP_GT: return a > b;
        case OP_GE: return a >= b;
        case OP_EQ: return a == b;
        case OP_NE: return a != b;
        case OP_AND: return a != 0 && b != 0;
        case OP_OR: return a != 0 || b != 0;
        default: return std::nan("");
    }
}

// ---- jets -------------------------------------------------------------------------------------------------------

// distance in x to the point where u reaches 0, to first order
inline double zero_dist(const Jet &u) {
    if (u.d1 != 0) return std::fabs(u.v / u.d1);
    return u.v == 0 ? 0 : INF;
}

// g(u) given g, g', g'' at u.v
inline Jet chain(const Jet &u, double g, double g1, double g2, double dist) {
    Jet r;
    r.v = g;
    r.d1 = g1 * u.d1;
    r.d2 = g2 * u.d1 * u.d1 + g1 * u.d2;
    r.dist = std::min(u.dist, dist);
    r.brk = u.brk;
    return r;
}

inline Jet jet_const(double c) {
    Jet r;
    r.v = c;
    return r;
}

inline Jet jet_mul(const Jet &a, const Jet &b) {
    Jet r;
    r.v = a.v * b.v;
    r.d1 = a.d1 * b.v + a.v * b.d1;
    r.d2 = a.d2 * b.v + 2 * a.d1 * b.d1 + a.v * b.d2;
    r.dist = std::min(a.dist, b.dist);
    r.brk = std::min(a.brk, b.brk);
    return r;
}

// distance in x to the point where a and b cross, to first order
inline double cross_dist(const Jet &a, const Jet &b) {
    Jet d;
    d.v = a.v - b.v;
    d.d1 = a.d1 - b.d1;
    return zero_dist(d);
}

// a condition or a choice between a and b: the result of the comparison, and where it flips
inline Jet jet_compare(OpCode op, const Jet &a, const Jet &b) {
    if (op == OP_MIN || op == OP_MAX) {
        Jet r = (op == OP_MIN) == (a.v < b.v) ? a : b;
        r.brk = std::min({a.brk, b.brk, cross_dist(a, b)});
        return r;
    }
    Jet r;
    r.v = binary_value(op, a.v, b.v);
    r.dist = std::min(a.dist, b.dist);
    r.brk = std::min(a.brk, b.brk);
    if (op != OP_AND && op != OP_OR) r.brk = std::min(r.brk, cross_dist(a, b));
    return r;
}

inline Jet jet_recip(const Jet &u) {
    double g = 1.0 / u.v;
    return chain(u, g, -g * g, 2 * g * g * g, zero_dist(u));
}

inline Jet jet_powc(const Jet &u, double p) {
    if (p == 2) return chain(u, u.v * u.v, 2 * u.v, 2, INF);
    bool smooth = p >= 0 && p == std::floor(p);
    double g = std::pow(u.v, p);
    double g1 = p * std::pow(u.v, p - 1);
    double g2 = p * (p - 1) * std::pow(u.v, p - 2);
    return chain(u, g, g1, g2, smooth ? INF : zero_dist(u));
}

inline Jet jet_unary(OpCode op, const Jet &u) {
    switch (op) {
        case OP_NEG: return chain(u, -u.v, -1, 0, INF);
        case OP_SIN: return chain(u, std::sin(u.v), std::cos(u.v), -std::sin(u.v), INF);
        case OP_COS: return chain(u, std::cos(u.v), -std::sin(u.v), -std::cos(u.v), INF);
        case OP_TAN: {
            double t = std::tan(u.v), s = 1 + t * t;
            double dist = u.d1 != 0 ? std::fabs(std::cos(u.v) / u.d1) : INF;
            return chain(u, t, s, 2 * t * s, dist);
        }
        case OP_ASIN: {
            double w = 1 - u.v * u.v, s = std::sqrt(w);
            double dist = u.d1 != 0 ? std::fabs((1 - std::fabs(u.v)) / u.d1) : INF;
            return chain(u, std::asin(u.v), 1 / s, u.v / (w * s), dist);
        }
        case OP_ACOS: {
            double w = 1 - u.v * u.v, s = std::sqrt(w);
            double dist = u.d1 != 0 ? std::fabs((1 - std::fabs(u.v)) / u.d1) : INF;
            return chain(u, std::acos(u.v), -1 / s, -u.v / (w * s), dist);
        }
        case OP_ATAN: {
            double w = 1 + u.v * u.v;
            return chain(u, std::atan(u.v), 1 / w, -2 * u.v / (w * w), INF);
        }
        case OP_SINH: return chain(u, std::sinh(u.v), std::cosh(u.v), std::sinh(u.v), INF);
        case OP_COSH: return chain(u, std::cosh(u.v), std::sinh(u.v), std::cosh(u.v), INF);
        case OP_TANH: {
            double t = std::tanh(u.v), s = 1 - t * t;
            return chain(u, t, s, -2 * t * s, INF);
        }
        case OP_EXP: {
            double g = std::exp(u.v);
            return chain(u, g, g, g, INF);
        }
        case OP_LOG: return chain(u, std::log(u.v), 1 / u.v, -1 / (u.v * u.v), zero_dist(u));
        case OP_LOG2: return chain(u, std::log2(u.v), 1 / (u.v * M_LN2), -1 / (u.v * u.v * M_LN2), zero_dist(u));
        case OP_LOG10: return chain(u, std::log10(u.v), 1 / (u.v * M_LN10), -1 / (u.v * u.v * M_LN10), zero_dist(u));
        case OP_SQRT: {
            double g = std::sqrt(u.v);
            return chain(u, g, 0.5 / g, -0.25 / (g * u.v), zero_dist(u));
        }
        case OP_CBRT: {
            double g = std::cbrt(u.v);
            return chain(u, g, g / (3 * u.v), -2 * g / (9 * u.v * u.v), zero_dist(u));
        }
        case OP_ABS: {
            double s = u.v < 0 ? -1 : 1;
            return chain(u, std::fabs(u.v), s, 0, INF);
        }
        case OP_SIGMOID: {
            double g = sigmoid(u.v), s = g * (1 - g);
            return chain(u, g, s, s * (1 - 2 * g), INF);
        }
        case OP_NOT: return chain(u, u.v == 0, 0, 0, INF);
        case OP_FLOOR:
        case OP_CEIL:
        case OP_ROUND:
        case OP_SIGN: {
            // constant around u, up to the nearest value of u where the step changes
            double g = unary_value(op, u.v), to_step;
            if (op == OP_SIGN) to_step = std::fabs(u.v);
            else if (op == OP_ROUND) to_step = 0.5 - std::fabs(u.v - g);
            else to_step = std::min(u.v - std::floor(u.v), std::ceil(u.v) - u.v);
            Jet r = chain(u, g, 0, 0, INF);
            double d = u.d1 != 0 ? std::fabs(to_step / u.d1) : (to_step == 0 ? 0 : INF);
            r.brk = std::min(r.brk, std::isfinite(u.v) ? d : INF);
            return r;
        }
        default: return jet_const(std::nan(""));
    }
}

// ---- syntax tree ------------------------------------------------------------------------------------------------

struct Node {
    OpCode op;
    double c = 0;
    std::unique_ptr<Node> a, b, c3;
};
using NodePtr = std::unique_ptr<Node>;

NodePtr make_const(double c) {
    auto n = std::make_unique<Node>();
    n->op = OP_CONST;
    n->c = c;
    return n;
}

NodePtr make_unary(OpCode op, NodePtr a) {
    if (a->op == OP_CONST) return make_const(unary_value(op, a->c));
    auto n = std::make_unique<Node>();
    n->op = op;
    n->a = std::move(a);
    return n;
}

NodePtr make_binary(OpCode op, NodePtr a, NodePtr b) {
    if (a->op == OP_CONST && b->op == OP_CONST) return make_const(binary_value(op, a->c, b->c));
    auto n = std::make_unique<Node>();
    n->op = op;
    n->a = std::move(a);
    n->b = std::move(b);
    return n;
}

NodePtr make_select(NodePtr cond, NodePtr a, NodePtr b) {
    if (cond->op == OP_CONST) return cond->c != 0 ? std::move(a) : std::move(b);
    auto n = std::make_unique<Node>();
    n->op = OP_SELECT;
    n->a = std::move(cond);
    n->b = std::move(a);
    n->c3 = std::move(b);
    return n;
}

class Parser {
   public:
    explicit Parser(const std::string &s) : s_(s) {}

    NodePtr parse() {
        NodePtr n = expr();
        skip();
        if (pos_ != s_.size()) fail("unexpected '" + std::string(1, s_[pos_]) + "'");
        return n;
    }

   private:
    const std::string &s_;
    size_t pos_ = 0;

    [[noreturn]] void fail(const std::string &what) const {
        throw std::invalid_argument("QPET: cannot parse QoI \"" + s_ + "\": " + what + " at position " +
                                    std::to_string(pos_));
    }

    void skip() {
        while (pos_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[pos_]))) pos_++;
    }

    char peek() {
        skip();
        return pos_ < s_.size() ? s_[pos_] : '\0';
    }

    bool accept(char c) {
        if (peek() != c) return false;
        pos_++;
        return true;
    }

    void expect(char c) {
        if (!accept(c)) fail(std::string("expected '") + c + "'");
    }

    // is the next token the word (such as "and"), not followed by more letters?
    bool at_word(const char *word) {
        skip();
        size_t n = std::strlen(word);
        if (pos_ + n > s_.size()) return false;
        for (size_t i = 0; i < n; i++)
            if (std::tolower(static_cast<unsigned char>(s_[pos_ + i])) != word[i]) return false;
        return !(pos_ + n < s_.size() && (std::isalnum(static_cast<unsigned char>(s_[pos_ + n])) || s_[pos_ + n] == '_'));
    }

    bool accept_word(const char *word) {
        if (!at_word(word)) return false;
        pos_ += std::strlen(word);
        return true;
    }

    // condition ? value : value, the lowest precedence; chains to the right
    NodePtr expr() {
        NodePtr n = disjunction();
        if (!accept('?')) return n;
        NodePtr a = expr();
        expect(':');
        return make_select(std::move(n), std::move(a), expr());
    }

    NodePtr disjunction() {
        NodePtr n = conjunction();
        while (accept_word("or")) n = make_binary(OP_OR, std::move(n), conjunction());
        return n;
    }

    NodePtr conjunction() {
        NodePtr n = negation();
        while (accept_word("and")) n = make_binary(OP_AND, std::move(n), negation());
        return n;
    }

    NodePtr negation() {
        if (accept_word("not")) return make_unary(OP_NOT, negation());
        return comparison();
    }

    NodePtr comparison() {
        NodePtr n = sum();
        char c = peek(), next = pos_ + 1 < s_.size() ? s_[pos_ + 1] : '\0';
        OpCode op;
        size_t len = 1;
        if (c == '<') op = next == '=' ? (len = 2, OP_LE) : OP_LT;
        else if (c == '>') op = next == '=' ? (len = 2, OP_GE) : OP_GT;
        else if (c == '=' && next == '=') op = (len = 2, OP_EQ);
        else if (c == '!' && next == '=') op = (len = 2, OP_NE);
        else return n;
        pos_ += len;
        return make_binary(op, std::move(n), sum());
    }

    NodePtr sum() {
        NodePtr n = term();
        for (;;) {
            if (accept('+')) n = make_binary(OP_ADD, std::move(n), term());
            else if (accept('-')) n = make_binary(OP_SUB, std::move(n), term());
            else return n;
        }
    }

    NodePtr term() {
        NodePtr n = unary();
        for (;;) {
            if (peek() == '*' && !(pos_ + 1 < s_.size() && s_[pos_ + 1] == '*')) {
                pos_++;
                n = make_binary(OP_MUL, std::move(n), unary());
            } else if (accept('/')) {
                n = make_binary(OP_DIV, std::move(n), unary());
            } else {
                return n;
            }
        }
    }

    NodePtr unary() {
        if (accept('-')) return make_unary(OP_NEG, unary());
        if (accept('+')) return unary();
        return power();
    }

    NodePtr power() {
        NodePtr base = primary();
        if (accept('^')) return make_binary(OP_POW, std::move(base), unary());
        if (peek() == '*' && pos_ + 1 < s_.size() && s_[pos_ + 1] == '*') {
            pos_ += 2;
            return make_binary(OP_POW, std::move(base), unary());
        }
        return base;
    }

    NodePtr primary() {
        char c = peek();
        if (c == '(') {
            pos_++;
            NodePtr n = expr();
            expect(')');
            return n;
        }
        if (c == '|') {
            pos_++;
            NodePtr n = expr();
            expect('|');
            return make_unary(OP_ABS, std::move(n));
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            const char *begin = s_.c_str() + pos_;
            char *end = nullptr;
            double v = std::strtod(begin, &end);
            if (end == begin) fail("bad number");
            pos_ += end - begin;
            NodePtr n = make_const(v);
            // implicit multiplication: 10x, 2(x+1), 3sin(x)
            char next = peek();
            if ((std::isalpha(static_cast<unsigned char>(next)) || next == '_' || next == '(') && !at_word("and") &&
                !at_word("or"))
                return make_binary(OP_MUL, std::move(n), power());
            return n;
        }
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            std::string id;
            while (pos_ < s_.size() && (std::isalnum(static_cast<unsigned char>(s_[pos_])) || s_[pos_] == '_'))
                id += static_cast<char>(std::tolower(static_cast<unsigned char>(s_[pos_++])));
            if (peek() == '(') return call(id);
            if (id == "x") {
                auto n = std::make_unique<Node>();
                n->op = OP_VAR;
                return n;
            }
            if (id == "pi") return make_const(M_PI);
            if (id == "e") return make_const(M_E);
            fail("unknown name '" + id + "'");
        }
        fail(c ? "unexpected '" + std::string(1, c) + "'" : "unexpected end");
    }

    NodePtr call(const std::string &id) {
        expect('(');
        NodePtr a = expr();
        NodePtr b, c;
        if (accept(',')) b = expr();
        if (b && accept(',')) c = expr();
        expect(')');
        if (id == "if") {  // if(condition, value if true, value if false)
            if (!c) fail("if takes three arguments");
            return make_select(std::move(a), std::move(b), std::move(c));
        }
        if (c) fail(id + " does not take three arguments");
        if (id == "min" || id == "max") {
            if (!b) fail(id + " takes two arguments");
            return make_binary(id == "min" ? OP_MIN : OP_MAX, std::move(a), std::move(b));
        }
        if (id == "pow") {
            if (!b) fail("pow takes two arguments");
            return make_binary(OP_POW, std::move(a), std::move(b));
        }
        if (id == "log" && b)  // log(x, base)
            return make_binary(OP_DIV, make_unary(OP_LOG, std::move(a)), make_unary(OP_LOG, std::move(b)));
        if (b) fail(id + " takes one argument");
        for (const Func &f : FUNCS)
            if (id == f.name) return make_unary(f.op, std::move(a));
        fail("unknown function '" + id + "'");
    }
};

// polynomial degree in x if it is 0 or 1, otherwise 2
int degree(const Node &n) {
    switch (n.op) {
        case OP_CONST: return 0;
        case OP_VAR: return 1;
        case OP_NEG: return degree(*n.a);
        case OP_ADD:
        case OP_SUB: return std::max(degree(*n.a), degree(*n.b));
        case OP_MUL: return std::min(2, degree(*n.a) + degree(*n.b));
        case OP_DIV: return n.b->op == OP_CONST ? degree(*n.a) : 2;
        default: return 2;
    }
}

}  // namespace

struct ExprCompiler {
    std::vector<Expr::Op> prog;
    size_t depth = 0, max_depth = 0;

    void push(OpCode op, double c, int stack_change) {
        prog.push_back({op, c});
        depth += stack_change;
        max_depth = std::max(max_depth, depth);
    }

    void emit(const Node &n) {
        switch (n.op) {
            case OP_CONST: push(OP_CONST, n.c, 1); return;
            case OP_VAR: push(OP_VAR, 0, 1); return;
            case OP_ADD:
            case OP_SUB:
            case OP_MUL:
            case OP_DIV:
            case OP_POW: emit_binary(n); return;
            case OP_MIN:
            case OP_MAX:
            case OP_LT:
            case OP_LE:
            case OP_GT:
            case OP_GE:
            case OP_EQ:
            case OP_NE:
            case OP_AND:
            case OP_OR:
                emit(*n.a);
                emit(*n.b);
                push(n.op, 0, -1);
                return;
            case OP_SELECT:
                emit(*n.a);
                emit(*n.b);
                emit(*n.c3);
                push(OP_SELECT, 0, -2);
                return;
            default:
                emit(*n.a);
                push(n.op, 0, 0);
                return;
        }
    }

    void emit_binary(const Node &n) {
        const Node &a = *n.a, &b = *n.b;
        if (b.op == OP_CONST) {
            emit(a);
            switch (n.op) {
                case OP_ADD: push(OP_ADDC, b.c, 0); return;
                case OP_SUB: push(OP_ADDC, -b.c, 0); return;
                case OP_MUL: push(OP_MULC, b.c, 0); return;
                case OP_DIV: push(OP_DIVC, b.c, 0); return;
                default: push(OP_POWC, b.c, 0); return;
            }
        }
        if (a.op == OP_CONST && (n.op != OP_POW || a.c > 0)) {
            emit(b);
            switch (n.op) {
                case OP_ADD: push(OP_ADDC, a.c, 0); return;
                case OP_SUB:
                    push(OP_NEG, 0, 0);
                    push(OP_ADDC, a.c, 0);
                    return;
                case OP_MUL: push(OP_MULC, a.c, 0); return;
                case OP_DIV: push(OP_RDIVC, a.c, 0); return;
                default:  // c^u = exp(u ln c)
                    push(OP_MULC, std::log(a.c), 0);
                    push(OP_EXP, 0, 0);
                    return;
            }
        }
        emit(a);
        emit(b);
        push(n.op, 0, -1);
    }
};

Expr::Expr(const std::string &text) : text_(text) {
    NodePtr root = Parser(text).parse();
    linear_ = degree(*root) <= 1;
    ExprCompiler c;
    c.emit(*root);
    if (c.max_depth > MAX_DEPTH) throw std::invalid_argument("QPET: QoI expression is nested too deeply");
    prog_ = std::move(c.prog);
    depth_ = c.max_depth;
}

namespace {

#define QPET_LOOP(expr) \
    for (size_t i = 0; i < n; i++) a[i] = (expr); \
    break

// st holds depth * n doubles
template <class Ops>
void run(const Ops &prog, const double *x, double *out, size_t n, double *st) {
    size_t sp = 0;
    for (const auto &op : prog) {
        double *a = st + (sp ? sp - 1 : 0) * n;
        const double c = op.c;
        switch (op.code) {
            case OP_CONST:
                a = st + sp++ * n;
                std::fill(a, a + n, c);
                break;
            case OP_VAR:
                a = st + sp++ * n;
                std::copy(x, x + n, a);
                break;
            case OP_ADD:
            case OP_SUB:
            case OP_MUL:
            case OP_DIV:
            case OP_POW:
            case OP_MIN:
            case OP_MAX:
            case OP_LT:
            case OP_LE:
            case OP_GT:
            case OP_GE:
            case OP_EQ:
            case OP_NE:
            case OP_AND:
            case OP_OR: {
                const double *b = a;
                a -= n;
                sp--;
                switch (op.code) {
                    case OP_ADD: QPET_LOOP(a[i] + b[i]);
                    case OP_SUB: QPET_LOOP(a[i] - b[i]);
                    case OP_MUL: QPET_LOOP(a[i] * b[i]);
                    case OP_DIV: QPET_LOOP(a[i] / b[i]);
                    case OP_POW: QPET_LOOP(std::pow(a[i], b[i]));
                    case OP_MIN: QPET_LOOP(a[i] < b[i] ? a[i] : b[i]);
                    case OP_MAX: QPET_LOOP(a[i] > b[i] ? a[i] : b[i]);
                    case OP_LT: QPET_LOOP(a[i] < b[i]);
                    case OP_LE: QPET_LOOP(a[i] <= b[i]);
                    case OP_GT: QPET_LOOP(a[i] > b[i]);
                    case OP_GE: QPET_LOOP(a[i] >= b[i]);
                    case OP_EQ: QPET_LOOP(a[i] == b[i]);
                    case OP_NE: QPET_LOOP(a[i] != b[i]);
                    case OP_AND: QPET_LOOP(a[i] != 0 && b[i] != 0);
                    default: QPET_LOOP(a[i] != 0 || b[i] != 0);
                }
                break;
            }
            case OP_SELECT: {  // both branches were computed for every value; keep the one the condition picks
                const double *no = a, *yes = a - n;
                a -= 2 * n;
                sp -= 2;
                QPET_LOOP(a[i] != 0 ? yes[i] : no[i]);
            }
            case OP_NOT: QPET_LOOP(a[i] == 0);
            case OP_FLOOR: QPET_LOOP(std::floor(a[i]));
            case OP_CEIL: QPET_LOOP(std::ceil(a[i]));
            case OP_ROUND: QPET_LOOP(std::round(a[i]));
            case OP_SIGN: QPET_LOOP(sign(a[i]));
            case OP_ADDC: QPET_LOOP(a[i] + c);
            case OP_MULC: QPET_LOOP(a[i] * c);
            case OP_DIVC: QPET_LOOP(a[i] / c);
            case OP_RDIVC: QPET_LOOP(c / a[i]);
            case OP_POWC:
                if (c == 2) {
                    QPET_LOOP(a[i] * a[i]);
                } else if (c == 3) {
                    QPET_LOOP(a[i] * a[i] * a[i]);
                } else {
                    QPET_LOOP(std::pow(a[i], c));
                }
            case OP_NEG: QPET_LOOP(-a[i]);
            case OP_SIN: QPET_LOOP(std::sin(a[i]));
            case OP_COS: QPET_LOOP(std::cos(a[i]));
            case OP_TAN: QPET_LOOP(std::tan(a[i]));
            case OP_ASIN: QPET_LOOP(std::asin(a[i]));
            case OP_ACOS: QPET_LOOP(std::acos(a[i]));
            case OP_ATAN: QPET_LOOP(std::atan(a[i]));
            case OP_SINH: QPET_LOOP(std::sinh(a[i]));
            case OP_COSH: QPET_LOOP(std::cosh(a[i]));
            case OP_TANH: QPET_LOOP(std::tanh(a[i]));
            case OP_EXP: QPET_LOOP(std::exp(a[i]));
            case OP_LOG: QPET_LOOP(std::log(a[i]));
            case OP_LOG2: QPET_LOOP(std::log2(a[i]));
            case OP_LOG10: QPET_LOOP(std::log10(a[i]));
            case OP_SQRT: QPET_LOOP(std::sqrt(a[i]));
            case OP_CBRT: QPET_LOOP(std::cbrt(a[i]));
            case OP_ABS: QPET_LOOP(std::fabs(a[i]));
            case OP_SIGMOID: QPET_LOOP(sigmoid(a[i]));
        }
    }
    std::copy(st, st + n, out);
}

#undef QPET_LOOP

}  // namespace

double Expr::eval(double x) const {
    double st[MAX_DEPTH], out;
    run(prog_, &x, &out, 1, st);
    return out;
}

void Expr::eval(const double *x, double *out, size_t n) const {
    thread_local std::vector<double> stack;
    if (stack.size() < depth_ * CHUNK) stack.resize(depth_ * CHUNK);
    for (size_t i = 0; i < n; i += CHUNK) {
        size_t m = std::min(CHUNK, n - i);
        run(prog_, x + i, out + i, m, stack.data());
    }
}

Jet Expr::jet(double x) const {
    Jet st[MAX_DEPTH];
    size_t sp = 0;
    for (const Op &op : prog_) {
        switch (op.code) {
            case OP_CONST: st[sp++] = jet_const(op.c); break;
            case OP_VAR:
                st[sp] = Jet();
                st[sp].v = x;
                st[sp].d1 = 1;
                sp++;
                break;
            case OP_ADD:
            case OP_SUB: {
                Jet &a = st[sp - 2];
                const Jet &b = st[sp - 1];
                double s = op.code == OP_ADD ? 1 : -1;
                a.v += s * b.v;
                a.d1 += s * b.d1;
                a.d2 += s * b.d2;
                a.dist = std::min(a.dist, b.dist);
                a.brk = std::min(a.brk, b.brk);
                sp--;
                break;
            }
            case OP_MIN:
            case OP_MAX:
            case OP_LT:
            case OP_LE:
            case OP_GT:
            case OP_GE:
            case OP_EQ:
            case OP_NE:
            case OP_AND:
            case OP_OR:
                st[sp - 2] = jet_compare(static_cast<OpCode>(op.code), st[sp - 2], st[sp - 1]);
                sp--;
                break;
            case OP_SELECT: {
                const double flips = st[sp - 3].brk;
                st[sp - 3] = st[sp - 3].v != 0 ? st[sp - 2] : st[sp - 1];
                st[sp - 3].brk = std::min(st[sp - 3].brk, flips);
                sp -= 2;
                break;
            }
            case OP_MUL:
                st[sp - 2] = jet_mul(st[sp - 2], st[sp - 1]);
                sp--;
                break;
            case OP_DIV:
                st[sp - 2] = jet_mul(st[sp - 2], jet_recip(st[sp - 1]));
                sp--;
                break;
            case OP_POW: {  // u^w = exp(w log u)
                Jet e = jet_mul(st[sp - 1], jet_unary(OP_LOG, st[sp - 2]));
                st[sp - 2] = jet_unary(OP_EXP, e);
                sp--;
                break;
            }
            case OP_ADDC: st[sp - 1].v += op.c; break;
            case OP_MULC:
                st[sp - 1].v *= op.c;
                st[sp - 1].d1 *= op.c;
                st[sp - 1].d2 *= op.c;
                break;
            case OP_DIVC:
                st[sp - 1].v /= op.c;
                st[sp - 1].d1 /= op.c;
                st[sp - 1].d2 /= op.c;
                break;
            case OP_RDIVC: {
                Jet r = jet_recip(st[sp - 1]);
                r.v *= op.c;
                r.d1 *= op.c;
                r.d2 *= op.c;
                st[sp - 1] = r;
                break;
            }
            case OP_POWC: st[sp - 1] = jet_powc(st[sp - 1], op.c); break;
            default: st[sp - 1] = jet_unary(static_cast<OpCode>(op.code), st[sp - 1]); break;
        }
    }
    return st[0];
}

double pointwise_eb(const Expr &f, double x, double tau, double eb_cap) {
    Jet j = f.jet(x);
    double a = std::fabs(j.d1), b = std::fabs(j.d2);
    if (!std::isfinite(j.v) || !std::isfinite(a) || !std::isfinite(b)) return 0;
    // largest e with a e + b e^2 / 2 <= tau, written so that it does not cancel when b is small
    double s = a + std::sqrt(a * a + 2 * b * tau);
    double eb = s > 0 ? 2 * tau / s : INF;
    eb = std::min({eb, j.dist, eb_cap});
    if (!(eb > 0)) return 0;
    auto error_at = [&](double e) { return std::max(std::fabs(f.eval(x - e) - j.v), std::fabs(f.eval(x + e) - j.v)); };
    // A piecewise function can jump where a condition flips. If stepping just across the nearest such point is out
    // of tolerance, the value has to stay on its side of it.
    // Beyond it another piece takes over, which the derivatives at x say nothing about: walk outwards from that
    // point in doubling steps and stop before the first one that is out of tolerance.
    if (j.brk < eb) {
        // "just across" and "just short of" leave room for the rounding of x +- e
        const double margin = std::max(j.brk, std::fabs(x)) * 1e-9 + std::numeric_limits<double>::min();
        const double across = j.brk + margin;
        if (!(error_at(across) <= tau)) {
            eb = j.brk - margin;
        } else {
            for (double e = 2 * across; e < eb; e *= 2) {
                if (!(error_at(e) <= tau)) {
                    eb = e / 2;
                    break;
                }
            }
        }
        if (!(eb > 0)) return 0;
    }
    if (!std::isfinite(eb)) return eb;
    // The Taylor estimate ignores the third derivative. Check both ends and shrink while one is out of tolerance:
    // in proportion to the excess when it is small, by half otherwise.
    for (int k = 0; k < 8; k++) {
        double err = error_at(eb);
        if (err <= tau) break;
        eb *= std::isfinite(err) ? std::max(0.5, 0.999 * tau / err) : 0.5;
    }
    return eb;
}

}  // namespace qpet
