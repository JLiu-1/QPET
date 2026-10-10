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
| step functions | `floor`, `ceil`, `round` (halves away from zero), `sign` (-1, 0 or 1) |
| functions of two arguments | `log(x, base)`, `pow(a, b)`, `min(a, b)`, `max(a, b)` |
| comparisons | `<`, `<=`, `>`, `>=`, `==`, `!=` |
| logic | `and`, `or`, `not` |
| conditional | `condition ? a : b`, or `if(condition, a, b)` |

Names are not case-sensitive. Spaces are ignored.

Precedence, from tightest: functions and parentheses; power (right-associative, so `2^-x` and `x^2^3 = x^(2^3)`);
unary minus (`-x^2` is `-(x^2)`); `*` and `/`; `+` and `-`. The implicit product binds as tightly as a function call:
`1/2x` is `1/(2x)`, and `2x^2` is `2(x^2)`.

`log` and `ln` are the natural logarithm. `sigmoid(x)` is `1/(1+exp(-x))`.

## Piecewise functions

A conditional picks one of two expressions by a condition on the value:

```
x < 0 ? 0 : sqrt(x)                  0 below zero, the square root above
x < 1 ? x^2 : x < 5 ? x : 5          three pieces: conditionals chain to the right
if(x >= 2 and x < 3, 10, x)          the same thing written as a function, with a combined condition
x > 273.15 ? 1 : 0                   a threshold: 1 where the value is above it
min(max(x, 0), 100)                  the value clipped to [0, 100]
```

The conditional has the lowest precedence, then `or`, `and`, `not`, the comparisons, and then arithmetic, so
`x < 1 ? x^2 : x + 1` needs no parentheses. Use parentheses to put a conditional inside a larger expression:
`2 * (x < 0 ? -x : x)`.

Only the piece that applies to a value has to be finite there: `x > 0 ? log2(x) : 0` is fine for values that are
not positive.

What the tolerance means at a breakpoint: it is the same guarantee as everywhere, `|f(x) - f(x')| <= tau`. If the
function jumps there by more than the tolerance, a value is never moved across the breakpoint, and a value exactly on
it is reproduced exactly. With `x > c ? 1 : 0` and a tolerance below 1, this keeps every value on its side of the
threshold `c`. If the pieces meet (or jump by less than the tolerance), values may cross.

## Step functions

`floor`, `ceil`, `round` and `sign` are constant between their breakpoints, so they are piecewise functions with many
pieces, and the same rule holds: a value is not moved across a step that is higher than the tolerance.

```
floor(x / 0.5)          the index of the bin of width 0.5 that the value falls in
round(x)                the nearest integer
sign(x)                 -1, 0 or 1
```

With a tolerance below 1, `floor(x / w)` keeps every value in its bin of width `w`: a histogram with those bins is
reproduced exactly. With a tolerance of 1 a value may move to a neighbouring bin, and so on. Values that sit close to
a bin edge tolerate little error, so the finer the bins, the more values are stored losslessly.

## Examples

```
x^2          x^3          x^2+x+1       sqrt(x)        1/x          (x+3)^-1
log2(x)      ln(x)        log(x, 10)    2^x            e^x          exp(-x)
sin(10x)     tanh(x)      |x|           sigmoid(x)     1/(1+exp(-x))
tanh(x/300)  log2(x+1)    pow(x, 2.5)   x^x
```

## What is and is not supported

- One variable. Functions of several fields are not supported.
- Smooth functions, `abs`, `min`, `max`, piecewise functions and step functions.
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
