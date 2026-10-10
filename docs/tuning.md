# Tuning

The backend takes one error bound for the whole array. A small one leaves few values that break the QoI tolerance
and have to be stored losslessly, but the backend compresses little; a large one does the opposite. The tuning looks
for the error bound at which the backend's stream plus the losslessly stored values are smallest.

It never affects the guarantees, only the compression ratio and the compression time.

## What it does

### 1. Point samples

About a million values are read at a regular stride over the whole input (`point_samples`). For each, the largest
error it tolerates is estimated:

1. From the first two derivatives `a = |f'(x)|`, `b = |f''(x)|`: `eb = 2 tau / (a + sqrt(a^2 + 2 b tau))`, the
   largest `e` with `a e + b e^2 / 2 <= tau`.
2. Not beyond the distance to the nearest singularity of `f` (for example 0 for `log`), nor beyond the data bound or
   the value range.
3. `f` is evaluated at `x - eb` and `x + eb`; while either end is out of tolerance, `eb` is reduced.
4. Where `f` or a derivative is not finite, the estimate is 0: the value has to be exact.

For a regional QoI over blocks of `n` values, `tau` is first multiplied by
`max(1, 2 sqrt(n / (2 ln(2 / (1 - 0.999)))))`: the tolerance of a single value when the errors of a block are
taken as independent (Theorem 5.4 of the paper with `c = 2` and confidence 0.999; about 4.1 for a 4x4x4 block).

The sorted estimates are the distribution that the candidates come from. Being spread evenly, the point samples
represent the whole input, however unevenly the sensitive values are distributed.

### 2. A sample to compress

**Blocks** (default). Contiguous blocks are copied out of the input.

- The block edge is `block_edge`. Each dimension is split into equal tiles of at most that edge (98 becomes 2 x 49),
  and a block is one tile.
- The number of blocks is `sample_rate` x input size / block size, kept between `min_blocks` and `max_blocks`.
- On a small input, `min_blocks` blocks of the full edge would be a large share of it. The edge is then reduced
  until the sample is about `sample_rate` of the input, but not below `min_edge`.
- The tiles are picked with a low-discrepancy (Halton) sequence, a different base per dimension, so that the blocks
  are spread over the input and every part of it is equally likely to be taken.
- For a regional QoI, block edges and origins are multiples of the region edge: the regions inside a block are
  regions of the input.
- If the blocks would amount to half the input or more, the whole input is used instead.

**Cells** (`independent_cell = c > 1`, the default for `zfp`). For a compressor that codes every cell of `c` values
per dimension on its own, single cells are taken, as many as `sample_rate` allows (tens of thousands), and laid side
by side in one array. They compress to exactly the bytes they take in the whole input, and they cover it far more
evenly than a few blocks can. With a regional QoI the cell edge becomes the least common multiple of `c` and the
region edge, and blocks are used instead if that exceeds 32.

### 3. The cost of a candidate

For one error bound, every sample block is compressed and decompressed with the backend, checked and corrected
exactly as the whole input will be, and the corrected values are encoded as they will be in the stream. The cost is

```
(backend bytes + weight * lossless bytes) / sample values
```

Two corrections make up for what blocks cannot see:

- `weight` is the share of values whose estimated error bound is below the candidate ("sensitive" values) in the
  whole input, divided by their share in the blocks. Sensitive values can be concentrated in a part of the input that
  the blocks over- or under-represent.
- If the blocks hold almost none of them although the input has them (their share in the input is more than 20
  times their share in the blocks), what the blocks measured says nothing. The lossless bytes are then taken as
  `rare_outlier_share` x (share of sensitive values in the input) x (bytes of a value + 1) per value.
- With `subtract_block_overhead`, the size of an empty block is subtracted from the backend bytes of each block.

### 4. The search

1. The first candidates are the values of the estimates' distribution at `quantiles`, without the zeros and the
   duplicates. If one is left, it is used without compressing anything.
2. Each is costed; the cheapest is the current best.
3. For `refine_rounds` rounds, the geometric midpoints between the best and its two neighbours are costed. The search
   moves to a midpoint that is cheaper, or else narrows the interval to the two midpoints. It stops early when the
   interval is narrower than 5 %.

With `pow2_error_bounds` (the default for `zfp`) the candidates are rounded down to powers of two, and instead of
step 3 the exponent is walked up, then down, one step at a time. A direction is abandoned after three steps in a row
that do not pay. The walk may go up to 4096 times the largest candidate, because a compressor that stays far inside
its tolerance can have its best setting above every estimate.

## Parameters

Each parameter is taken from the first place that sets it: the user, the backend, the general default.

| name | meaning | general default |
|---|---|---|
| `sample_rate` | fraction of the input compressed per candidate | 0.01 |
| `block_edge` | edge of a sample block | 262144 (1-D), 512 (2-D), 64 (3-D), 24 (4-D) |
| `min_edge` | smallest edge the blocks shrink to on a small input | `block_edge`, 128 in 2-D |
| `min_blocks`, `max_blocks` | bounds on the number of blocks | 5 and 64; 12 and 64 in 2-D |
| `independent_cell` | sample single cells of this edge instead of blocks; `1` turns it off | off |
| `pow2_error_bounds` | `1`: try only powers of two | 0 |
| `refine_rounds` | rounds of narrowing around the best candidate | 4 |
| `subtract_block_overhead` | `1`: subtract the size of an empty block from every block | 0 |
| `rare_outlier_share` | share of the sensitive values assumed stored losslessly when the blocks missed them | 0.1 |
| `quantiles` | first candidates, as quantiles of the estimates | 0, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1 |
| `point_samples` | values the estimates are computed for | 1048576 |

What the built-in backends set:

| backend | setting | why |
|---|---|---|
| `zfp` | `independent_cell=4`, `pow2_error_bounds=1`, `rare_outlier_share=0.01` | zfp codes each 4^d cell on its own, only reads the binary exponent of the tolerance, and its errors stay far below it |
| `sperr` | `rare_outlier_share=0.02` | most of its errors are a small fraction of the tolerance |
| `sz3`, `szo` | none | |

How to set them:

| interface | form |
|---|---|
| command line | `-S min_blocks=24 -S quantiles=0:0.01:0.1:1` |
| C++ | `conf.tuning.min_blocks = 24;` or `conf.set("tuning.min_blocks", "24");` |
| C | `tuning.min_blocks=24` in the options text |
| Python | `tuning={"min_blocks": 24, "quantiles": [0, 0.01, 0.1, 1]}` |

In text form the quantiles are separated by `:`.

## Skipping the tuning

- `backend_eb` (`-E`): run the backend with this absolute error bound.
- `quantile` (`-Q`): use this quantile of the estimates.

Either way the estimate pass still runs, and the checks and corrections are unchanged.

## Looking at what it did

`-v` (or `verbose`) prints a line per candidate:

```
QPET tuning: eb 0.0183366   backend 0.4412 + lossless 0.0498 x 0.93 = 0.4874 bytes/value
             (stored losslessly: 1.1% of the sample, sensitive: 12.4% of the sample, 11.5% of the input)
```

The report of the compression shows how the choice turned out on the whole input: `outliers`, `outlier bytes` and
`backend bytes`. If the outliers take most of the stream, the error bound was too large; if there are none and the
ratio is low, it may have been too small.

## What is known not to work well

Measured during development, by comparing the tuning with a scan of the error bound on the whole input (20 cases
per backend):

- **Small blocks cost a compressor more than the same data inside a large array**, by an amount that grows with the
  error bound. This pushes the choice towards smaller error bounds than the best one.
- **Rare, clustered sensitive values** (a few storms in a precipitation field) can be missed by all blocks. The
  `rare_outlier_share` fallback limits the damage but is a fixed guess.
- **Few blocks give a noisy estimate.** On 2-D inputs, 48 blocks instead of 12 (`-S min_blocks=48`) roughly halved the
  loss against the best error bound for `sz3` and `szo`, for about three times the tuning time.
- **Very high compression ratios** (a stream of a few kilobytes): the size is no longer a smooth function of the
  error bound and a sample cannot predict it.
