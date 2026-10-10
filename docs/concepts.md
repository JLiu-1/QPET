# Concepts

## The problem

An error-bounded lossy compressor bounds the error of the data: every decompressed value `x'` is within `eb` of the
original `x`. An analysis usually looks at something computed from the data, a quantity of interest (QoI), and a bound
on the data says little about the error of a non-linear QoI. QPET compresses so that the QoI itself stays within a
tolerance, and does so on top of an existing compressor, without modifying it.

## What is guaranteed

Let `f` be the QoI function, a function of one value (see [QoI expressions](qoi.md)).

**Point-wise QoI** (the default). For every value:

```
|f(x) - f(x')| <= tau
```

**Regional QoI** (`region = r > 1`). The array is cut into blocks of `r` values per dimension, starting at the origin.
Blocks at the far ends are smaller when a dimension is not a multiple of `r`. For every block `B` with `n` values:

```
| (1/n) sum_{i in B} f(x_i)  -  (1/n) sum_{i in B} f(x'_i) |  <= tau
```

Individual values of `f` may be off by more than `tau`; only the block means are bound.

In both cases:

- **Data bound** (optional). If one is given, `|x - x'| <= eps` for every value.
- **Where `f` is not finite.** If `f(x)` is NaN or infinite (for example `log2` of a value that is not positive), the
  value is reproduced exactly: `x' == x`. In a regional QoI such a value adds nothing to the sums (`n` is still the
  number of values in the block). A NaN in the data is reproduced as a NaN.

These hold whatever the backend does. QPET decompresses what the backend produced, checks it against the original,
and stores every value that has to be exact losslessly. The guarantee therefore does not rely on the backend
respecting its error bound, nor on the estimates QPET makes along the way.

The errors are computed in double precision from the stored values. An analysis that evaluates the same QoI in single
precision, or with a different formula that is equal only mathematically, can see a difference of the order of the
rounding error on top of `tau`.

## Tolerances: absolute and relative

Both `tau` and `eps` can be given as absolute numbers or relative to a value range:

| | absolute (`ABS`) | relative (`REL`) |
|---|---|---|
| QoI tolerance, point-wise | `tau` | `tau * (max f(x) - min f(x))` over the input, taken where `f` is finite |
| QoI tolerance, regional | `tau` | `tau * (max - min)` of the block means of `f` over the input |
| data bound | `eps` | `eps * (max x - min x)` over the input |

If the range is zero (constant data), the relative value is used as an absolute one. The absolute values that were
used are reported by the compression and stored in the stream.

A relative tolerance depends on the extremes of the input. A few fill values such as `1e36` make the range, and with
it the tolerance, enormous: give an absolute tolerance for such data.

## How the compression works

1. **Estimate.** For a sample of the values, estimate the largest error each one tolerates: the second-order Taylor
   bound of the QPET paper (Theorems 5.1 and 5.2), then a check of `f` at both ends of the interval. For a regional
   QoI the tolerance of a single value is first enlarged by the factor of Theorem 5.4, because errors partly cancel
   in a mean.
2. **Tune.** The backend takes one error bound for the whole array. Candidates are taken from the distribution of
   the estimates. Each candidate is tried on a small sample of the input: compress, decompress, check, and count the
   bytes of the backend's stream plus those of the values that would be stored losslessly. The cheapest wins. See
   [Tuning](tuning.md).
3. **Compress** the whole input with the backend at that error bound.
4. **Check and correct.** Decompress, and compare with the original.
   - Point-wise: a value whose QoI error or data error is out of bounds is stored losslessly.
   - Regional: in a block whose mean is out of tolerance, values are restored one at a time, those that contribute
     most to the error of the mean first, until the block is within tolerance.
5. **Write** the stream: a header, the backend's stream, and the losslessly stored values ("outliers").

Decompression runs the backend's decompression and writes the stored values over the result. It does not need the
original data or the QoI.

Steps 1 and 2 only affect the compression ratio. If they choose badly, more values end up stored losslessly; the
bounds still hold.

## Cost

- Compression runs the backend's compression once and its decompression once on the whole input, plus the tuning on
  the sample, and reads the input a few more times for the ranges and the check. It needs memory for one decompressed
  copy of the input on top of what the backend uses.
- Decompression costs the backend's decompression plus the time to write the outliers.

## Limits

- One QoI function of one value. QoIs of several fields (for example a velocity magnitude from three components) are
  not supported.
- Regions are aligned blocks of one size. Sliding windows, masks and regions that span files are not supported.
- One error bound per array is given to the backend. The paper's SZ3-QPET and HPEZ-QPET also quantize a separate
  error bound for each value inside the compressor, which this layer cannot do; where the tolerance of the values
  varies strongly in space, they compress better.
- `float32` and `float64` arrays with 1 to 4 dimensions, held in memory.
