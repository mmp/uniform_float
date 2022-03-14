// uniformfloat.h
//
// Routines described in https://pharr.org/matt/blog/2022/03/14/sampling-float-intervals.html
//
// Releated previous work:
// - Christoph Conrads (2018): https://gitlab.com/christoph-conrads/rademacher-fpl
// - Olaf Bernstein (2021): https://github.com/camel-cdr/cauldron
//
// Licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef UNIFORMFLOAT_H
#define UNIFORMFLOAT_H

#include <cassert>
#include <cstdint>
#include <random>
#include <version>

#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

// Expect the user to provide implemetnations of these.
extern HD uint64_t Random64Bits();
extern HD uint32_t Random32Bits();

#ifdef __cpp_lib_bit_cast
  #include <bit>

  HD uint32_t ToBits(float f) { return std::bit_cast<uint32_t>(f); }
  HD float FromBits(uint32_t b) { return std::bit_cast<float>(b); }

#else
  #include <cstring>

  template <typename To, typename From> HD To bit_cast(From f) {
      static_assert(sizeof(To) == sizeof(From), "To and From types must have same size");
      To to;
      std::memcpy(&to, &f, sizeof(to));
      return to;
  }
  HD uint32_t ToBits(float f) { return bit_cast<uint32_t>(f); }
  HD float FromBits(uint32_t b) { return bit_cast<float>(b); }

#endif

namespace uniform_float {

// __clang__/__GNUC__ specializations via
// http://marc-b-reynolds.github.io/math/2020/02/16/LogOfUniform.html#helper-functions
// See further discussion there...
#if defined(__CUDA_ARCH__)
inline __device__ int CountLeadingZeros(uint64_t x) { return __clzll(x); }
#elif defined(__clang__)
inline int CountLeadingZeros(uint64_t x) { return x ? __builtin_clzl(x) : 64; }
#elif defined(__GNUC__)
inline int CountLeadingZeros(uint64_t x) { return __builtin_clzl(x); }
#else
#include <immintrin.h>
inline int CountLeadingZeros(uint64_t x) { return _lzcnt_u64(x); }
#endif

inline HD int SignBit(float f) { return ToBits(f) >> 31; }
inline HD int Exponent(float f) { return ((ToBits(f) >> 23) & 0xff) - 127; }
constexpr int SignificandMask = (1 << 23) - 1;
inline HD int Significand(float f) { return ToBits(f) & SignificandMask; }

inline HD uint32_t RandomSignificand() { return Random32Bits() & SignificandMask; }

inline HD float Float32FromParts(int sign, int exponent, int significand) {
    assert(sign == 0 || sign == 1);
    assert(exponent >= -127 && exponent <= 127);
    assert(significand >= 0 && significand < (1 << 23));
    return FromBits((sign << 31) | ((exponent + 127) << 23) | significand);
}

inline HD float FloatPow2(int exponent) {
    assert(exponent >= -126 && exponent <= 127);
    return FromBits((exponent + 127) << 23);
}

// Sample uniformly and comprehensively in [0,1)
// This is an implementation of
// http://marc-b-reynolds.github.io/distribution/2017/01/17/DenseFloat.html#the-parts-im-not-tell-you
inline HD float Sample01() {
    uint64_t bits = Random64Bits();
    int significand = bits & SignificandMask;
    if (int lz = CountLeadingZeros(bits); lz <= 40)
        return Float32FromParts(0, -1 - lz, significand);
    return 0x1p-64f * significand;
}

// Sample an exponent for [0,2^exponent)
inline HD int SampleToPowerOfTwoExponent(int exponent) {
    assert(exponent >= -127 && exponent <= 127);
    while (exponent > -126) {
        if (int lz = CountLeadingZeros(Random64Bits()); lz == 64)
            exponent -= 64;
        else
            return std::max(-127, exponent - 1 - lz);
    }
    return -127;
}

// Sample uniformly and comprehensively in [0,2^exponent)
inline HD float SampleToPowerOfTwo(int exponent) {
    int ex = SampleToPowerOfTwoExponent(exponent);
    return Float32FromParts(0, ex, RandomSignificand());
}

// Sample uniformly in [0,2^exponent) using a fixed number of bits
inline HD float SampleToPowerOfTwoFast(int exponent, uint64_t bits) {
    int significand = bits & SignificandMask;
    int lz = CountLeadingZeros(bits);
    if (lz == 41 && exponent - 41 > -127)
        return significand * 0x1p-23f * FloatPow2(exponent - 41);
    int ex = exponent - 1 - lz;
    return Float32FromParts(0, std::max(-127, ex), significand);
}

inline HD int SampleExponent(int emin, int emax) {
    int e = 0;
    while (true) {
        if (int lz = CountLeadingZeros(Random64Bits()); lz == 64)
            e += 64;
        else
            return emax - 1 - ((e + lz) % (emax - emin));
    }
}

// Sample uniformly and comprehensively in [2^emin, 2^emax).
inline HD float SampleExponentRange(int emin, int emax) {
    assert(emax > emin);
    int significand = RandomSignificand();
    return Float32FromParts(0, SampleExponent(emin, emax), significand);
}

// Sample uniformly in the range [a,b).
inline HD float SampleRange(float a, float b) {
    assert(a < b && a >= 0.f && b >= 0.f);
    int ea = Exponent(a), eb = Exponent(b);
    if (Significand(b) != 0) ++eb;
    while (true) {
       int e = (ea == -127) ? SampleToPowerOfTwoExponent(eb) :
                              SampleExponent(ea, eb);
       float v = Float32FromParts(0, e, RandomSignificand());
       if (v >= a && v < b)
           return v;
    }
}

}  // namespace uniform_float

#endif // UNIFORMFLOAT_H
