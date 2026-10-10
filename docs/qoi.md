# QoI expressions

The QoI is a function `f(x)` of one data value, written as text. It is compiled once into a small program; its
derivatives, which the estimates need, come from automatic differentiation, so there is nothing to derive by hand.

## Syntax

| element | forms |
|---|---|
| variable | `x` |
| numbers | `2`, `0.5`, `1e-3` |
| constants | `pi`, `e` |
| operators | `+`, `-`, `*`, `/`, `^` or `**` (power), unary `-` |
| grouping | `( )`, and `\|...\|` for the absolute value |
| implicit product | a number directly in front of a name or a parenthesis multiplies: `10x`, `2(x+1)`, `3sin(x)` |
| functions of one argument | `sin cos tan asin acos atan sinh cosh tanh exp log ln log2 log10 sqrt cbrt abs sigmoid` |
| functions of two arguments | `log(x, base)`, `pow(a, b)` |

Names are not case-sensitive. Spaces are ignored.

Precedence, from tightest: functions and parentheses; power (right-associative, so `2^-x` and `x^2^3 = x^(2^3)`);
unary minus (`-x^2` is `-(x^2)`); `*` and `/`; `+` and `-`. The implicit product binds as tightly as a function call:
`1/2x` is `1/(2x)`, and `2x^2` is `2(x^2)`.

`log` and `ln` are the natural logarithm. `sigmoid(x)` is `1/(1+exp(-x))`.

## Examples

```
x^2          x^3          x^2+x+1       sqrt(x)        1/x          (x+3)^-1
log2(x)      ln(x)        log(x, 10)    2^x            e^x          exp(-x)
sin(10x)     tanh(x)      |x|           sigmoid(x)     1/(1+exp(-x))
tanh(x/300)  log2(x+1)    pow(x, 2.5)   x^x
```

## What is and is not supported

- One variable. Functions of several fields are not supported.
- Smooth functions and `abs`. There are no comparisons, conditionals, `min`/`max`, rounding or thresholds.
- An expression that is nested more than 64 levels deep is rejected.

## Where the function is not finite

A value where `f` is NaN or infinite, such as `log2(x)` for `x <= 0` or `1/x` at 0, is reproduced exactly. For
`sqrt(x)` the function is finite at 0 but its slope is not; values at and around 0 simply end up with a very small
tolerance. If a large share of the data lies where `f` is not finite, all of it is stored losslessly and the
compression ratio suffers: consider shifting the argument (`log2(x+1)`), if that is what the analysis computes.

## How the expression is evaluated

In double precision, whatever the data type of the input, with the C math library. The same evaluation is used when
compressing and by `verify`. Two expressions that are equal mathematically, such as `tanh(x)` and
`(e^x-e^-x)/(e^x+e^-x)`, can differ in the last bits, and so can a tolerance check that sits exactly on the limit.

## Linear functions

For `f(x) = a x + b` the QoI error is `|a|` times the data error and every value tolerates the same error,
`min(tau/|a|, eps)`. With a point-wise QoI there is then a single candidate for the backend's error bound and no
tuning is needed (for `zfp`, whose errors stay far below its tolerance, the tuning still looks for a larger one).
