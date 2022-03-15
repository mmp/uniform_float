uniform_float
=============

This is a small header-only library implementing a handful of routines for
uniform sampling in floating-point.  See the blog posts [Sampling in
Floating Point (1/3): The Unit
Interval](https://pharr.org/matt/blog/2022/03/05/sampling-fp-unit-interval.html)
and [Sampling in Floating Point (2/3): 1D
Intervals](https://pharr.org/matt/blog/2022/03/14/sampling-float-intervals.html)
for further context and documentation.

The main entrypoints are:
- `Sample01()`: returns uniform floating-point values in [0,1).
- `SampleToPowerOfTwo(x)`: uniformly samples in [0,2^x).
- `SampleRange(a, b)`: uniformly samples in [a,b).

It is expected that the user will provide implementations of two functions
to generate random values:
```c++
uint64_t Random64Bits();
uint32_t Random32Bits();
```

For previous implementations of these ideas, see:
- Christoph Conrads's [Rademacher Floating Point
Library](https://gitlab.com/christoph-conrads/rademacher-fpl), in
particular the
[make_uniform_random_value()](https://gitlab.com/christoph-conrads/rademacher-fpl/-/blob/master/include/rademacher-fpl/impl/uniform-real-distribution.hpp#L225)
function.
- Olaf Bernstein's
[dist_uniformf_dense()](https://github.com/camel-cdr/cauldron/blob/7d5328441b1a1bc8143f627aebafe58b29531cb9/cauldron/random.h#L1604).
- Marc Reynolds's [Higher density uniform
  floats](http://marc-b-reynolds.github.io/distribution/2017/01/17/DenseFloat.html#the-parts-im-not-tell-you)
  blog post.
