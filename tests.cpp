
#include "uniformfloat.h"

#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <vector>

std::minstd_rand gen(6502);

uint64_t Random64Bits() {
    return gen();
}

uint32_t Random32Bits() {
    return gen();
}

#define CHECK(x) \
    do { if (!(x)) { fprintf(stderr, "  CHECK: " #x " failed. (Line %d).\n", __LINE__); ++nFailures; } } while (false)

#if 0
inline float SampleSameExponent(float a, float b) {
    assert(a < b && Exponent(a) == Exponent(b));
    int sa = Significand(a), sb = Significand(b);
    int sig = sa + RandomInt(sb - sa - 1);
    return Float32FromParts(SignBit(a), Exponent(a), sig);
}
#endif

using namespace uniform_float;

template <typename S>
int CheckPow2Sampler(S sampler, float min, float max) {
    // 1. Make sure exponents are roughly geometrically-distributed
    std::vector<int> exponentCount(255, 0);
    // 2. Check for uniformity by splitting [0,1) into 1024 buckets and
    // make sure each one has about the same number of entries.
    constexpr int nBuckets = 1024;
    std::vector<int> bucketCount(nBuckets, 0);

    assert(max > 0 && Significand(max) == 0);
    int maxExponent = Exponent(max) - 1;

    assert(gen.min() == 0);
    assert(gen.max() == 0xffffffffffffffff);

    constexpr int nUniform = 256*1024*1024;
    for (int i = 0; i < nUniform; ++i) {
        float u = sampler();
        assert(u >= min && u < max);  // assert intentional; we can't continue
        ++bucketCount[std::min<int>((u - min) / (max - min) * nBuckets, nBuckets - 1)];

        int e = Exponent(u);
        assert(e <= maxExponent);
        ++exponentCount[-e + maxExponent];
    }

    int nFailures = 0;
    if (min == 0) {
        float frac = 0.5f;
        for (int i = 0; i < exponentCount.size(); ++i, frac *= 0.5f) {
            //fprintf(stderr, "%d\t%d\t%d\t%d\t%.9g\n", maxExponent, maxExponent-i, i, exponentCount[i], nUniform*frac);
            if (maxExponent - i == -127)
                frac *= 2.f;  // denorms should have the same fraction as the next interval up

            if (maxExponent - i < -127)   // not a valid float at this point
                CHECK(exponentCount[i] == 0);
            else if (i < 10)
                CHECK(exponentCount[i] > 0.975 * frac * nUniform &&
                      exponentCount[i] < 1.025 * frac * nUniform);
            else if (i < 21)
                CHECK(exponentCount[i] > 0.75 * frac * nUniform &&
                      exponentCount[i] < 1.25 * frac * nUniform);
            else if (i < 29)
                CHECK(std::abs(exponentCount[i] - int(frac * nUniform)) <= 20);
            else
                CHECK(exponentCount[i] <= 3);
        }
    } else {
    }

    for (int i = 0; i < nBuckets; ++i) {
        //printf("%d\t%d\t%d\n", i, bucketCount[i], nUniform / nBuckets);
        CHECK(bucketCount[i] >= .98 * nUniform / nBuckets &&
              bucketCount[i] <= 1.02 * nUniform / nBuckets);
    }

    return nFailures;
}

template <typename S>
int CheckGeneralSampler(S sampler, float min, float max) {
    // 1. Check for uniformity by splitting [0,1) into 1024 buckets and
    // make sure each one has about the same number of entries.
    constexpr int nBuckets = 1024;
    std::vector<int> bucketCount(nBuckets, 0);

    constexpr int nUniform = 256*1024*1024;
    for (int i = 0; i < nUniform; ++i) {
        float u = sampler();
        assert(u >= min && u < max);  // assert intentional; we can't continue
        ++bucketCount[std::min<int>((u - min) / (max - min) * nBuckets, nBuckets - 1)];
    }

    int nFailures = 0;

    for (int i = 0; i < nBuckets; ++i) {
        //printf("%d\t%d\t%d\n", i, bucketCount[i], nUniform / nBuckets);
        CHECK(bucketCount[i] >= .98 * nUniform / nBuckets &&
              bucketCount[i] <= 1.02 * nUniform / nBuckets);
    }

    return nFailures;
}
int main() {
    using namespace uniform_float;
    int nFailures = 0;

    fprintf(stderr, "Checking basics\n");
    CHECK(CountLeadingZeros(0b1000000000000000000000000000000000000000000000000000000000000000) == 0);
    CHECK(CountLeadingZeros(0b0100000000000000000000000000000000000000000000000000000000000000) == 1);
    CHECK(CountLeadingZeros(0b0110000000000000000000000000000000000000000000000000000000000000) == 1);
    CHECK(CountLeadingZeros(0b0000110000000000000000000000000000000000000000000000000000000000) == 4);
    CHECK(CountLeadingZeros(0b0000000000000000000000000000000000000000000000000000000000000010) == 62);
    CHECK(CountLeadingZeros(0b0000000000000000000000000000000000000000000000000000000000000001) == 63);
    CHECK(CountLeadingZeros(0) == 64);

    // Deconstruct and some floats
    CHECK(SignBit(1.f) == 0);
    CHECK(SignBit(-1.f) == 1);
    CHECK(SignBit(-1245.f) == 1);

    CHECK(Exponent(1.f) == 0);
    CHECK(Exponent(0.5f) == -1);
    CHECK(Exponent(0.75f) == -1);
    CHECK(Exponent(-1.f) == 0);
    CHECK(Exponent(-0.5f) == -1);
    CHECK(Exponent(-0.75f) == -1);
    CHECK(Exponent(0.f) == -127);

    float f = 1.f;
    CHECK(f == Float32FromParts(SignBit(f), Exponent(f), Significand(f)));
    f = -1.f;
    CHECK(f == Float32FromParts(SignBit(f), Exponent(f), Significand(f)));
    f = 3.14159f;
    CHECK(f == Float32FromParts(SignBit(f), Exponent(f), Significand(f)));
    f = -f;
    CHECK(f == Float32FromParts(SignBit(f), Exponent(f), Significand(f)));

    for (int ex = -126; ex <= 127; ++ex)
        CHECK(FloatPow2(ex) == std::pow(2.f, ex));

    CHECK(Float32FromParts(0, -127, 1) == std::pow(2.f, -149)); // smallest denorm
    CHECK(Float32FromParts(1, -127, 1) == -std::pow(2.f, -149));

    // the width of the denorms, [0,2^-126] should be the same as the width
    // of the first power-of-2 sized span of normal numbers
    CHECK(Float32FromParts(0, -126, 0) == Float32FromParts(0, -125, 0) - Float32FromParts(0, -126, 0));

    ///////////////////////////////////////////////////////////////////////////
    // Test specific sampling routines

    // Sample01
    fprintf(stderr, "Checking Sample01()\n");
    nFailures += CheckPow2Sampler([]() { return Sample01(); }, 0.f, 1.f);

    // SampleToPowerOfTwo
    fprintf(stderr, "Checking SampleToPowerOfTwo()\n");
    for (int exponent : { 0, -1, -10, 15, 70, -90, -124 }) {
        auto sample = [&]() { return SampleToPowerOfTwo(exponent); };
        nFailures += CheckPow2Sampler(sample, 0.f, std::pow(2.f, exponent));
    }

    // SampleToPowerOfTwoFast
    fprintf(stderr, "Checking SampleToPowerOfTwoFast()\n");
    for (int exponent : { 0, -1, -10, 15, 70, -90, -124 }) {
        auto sample = [&]() { return SampleToPowerOfTwoFast(exponent, gen()); };
        nFailures += CheckPow2Sampler(sample, 0.f, std::pow(2.f, exponent));
    }

    // SampleExponent
    fprintf(stderr, "Checking SampleExponent()\n");
    // [1,8) -> expect 1/7, 2/7, 3/7
    int counts[3] = {0};
    for (int i = 0; i < 256*1024*1024; ++i) {
        int ex = SampleExponent(0, 3);
        assert(ex >= 0 && ex < 3);
        ++counts[ex];
    }
    for (int i = 0; i < 3; ++i)
        CHECK(counts[i] >= .98 * (i + 1.) / 3. &&
              counts[i] >= 1.02 * (i + 1.) / 3.);

    // SampleExponentRange
    fprintf(stderr, "Checking SampleExponentRange()\n");
    {
        auto sample = [&]() { return SampleExponentRange(1, 4); };
        nFailures += CheckPow2Sampler(sample, 2.f, 16.f);
    }
    {
        int e0 = -126, e1 = -110;
        auto sample = [&]() { return SampleExponentRange(e0, e1); };
        nFailures += CheckPow2Sampler(sample, std::pow(2.f, e0), std::pow(2.f, e1));
    }
    {
        int e0 = -120, e1 = 120;
        auto sample = [&]() { return SampleExponentRange(e0, e1); };
        nFailures += CheckPow2Sampler(sample, std::pow(2.f, e0), std::pow(2.f, e1));
    }
    {
        int e0 = -5, e1 = 5;
        auto sample = [&]() { return SampleExponentRange(e0, e1); };
        nFailures += CheckPow2Sampler(sample, std::pow(2.f, e0), std::pow(2.f, e1));
    }
    {
        int e0 = 15, e1 = 17;
        auto sample = [&]() { return SampleExponentRange(e0, e1); };
        nFailures += CheckPow2Sampler(sample, std::pow(2.f, e0), std::pow(2.f, e1));
    }

    // SampleRange
    fprintf(stderr, "Checking SampleRange()\n");
    {
        auto sample = [&]() { return SampleRange(0.f, 1.f); };
        nFailures += CheckGeneralSampler(sample, 0.f, 1.f);
    }
    {
        float min = std::pow(2.f, -125);
        auto sample = [&]() { return SampleRange(min, 1.f); };
        nFailures += CheckGeneralSampler(sample, min, 1.f);
    }
    {
        float max = std::pow(2.f, -125);
        auto sample = [&]() { return SampleRange(0.f, max); };
        nFailures += CheckGeneralSampler(sample, 0.f, max);
    }
    {
        float min = 1.5f, max = 8.5f;
        auto sample = [&]() { return SampleRange(min, max); };
        nFailures += CheckGeneralSampler(sample, min, max);
    }
    {
        float min = 0.f, max = 263.125f;
        auto sample = [&]() { return SampleRange(min, max); };
        nFailures += CheckGeneralSampler(sample, min, max);
    }
    {
        float min = std::pow(2.f, -126), max = 1.4426221 * std::pow(2.f, -124);
        auto sample = [&]() { return SampleRange(min, max); };
        nFailures += CheckGeneralSampler(sample, min, max);
    }

    return nFailures > 0;
}
