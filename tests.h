#pragma once

#include <cmath>
#include <iostream>
#include <set>
#include <vector>
#include "lfsr_hash.h"

#include <chrono>

namespace timer_n {

class Timer {
    std::chrono::time_point<std::chrono::steady_clock> mBegin;
public:
    Timer(): mBegin(std::chrono::steady_clock::now()) {};
    double elapsed_ns() {
        auto mEnd = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(mEnd - mBegin).count();
    }
    void reset() {
        mBegin = std::chrono::steady_clock::now();
    }
};

}

template <size_t N>
double bench(timer_n::Timer& timer, lfsr_hash::gens& g, const std::vector<uint8_t>& data) {
    if (data.size() > N) {
        return 1./0.;
    }
    timer.reset();
    [[maybe_unused]]auto hash128 = lfsr_hash::hash128<N>(g, data.data(), {1, 2, 3});
    const auto dt = timer.elapsed_ns();
    const double perf = (1000. * N) / dt;
    return perf;
}

static inline void lfsr_hash_benchmark() {
    std::cout << "Wait for 128-bit LFSR hash benchmark...\n";
    std::cout << std::flush;
    timer_n::Timer timer;
    lfsr_hash::gens g;
    constexpr size_t N = 4*1024*1024;
    auto v = std::vector<uint8_t>(N);
    std::cout << "Input array of " << N << " bytes is allocated.\n";
    std::vector<double> perfs {};
    constexpr size_t M = 128;
    for (int i=0; i<M; ++i) {
        const auto perf = bench<N>(timer, g, v);
        if (!std::isnan(perf) && !std::isinf(perf)) {
            perfs.push_back( perf );
        }
    }
    std::sort(perfs.begin(), perfs.end());
    std::cout << "128-bit LFSR hash median performance: " << perfs[perfs.size()/2] << " MB/s\n";
    std::cout << " All Ok! Completed.\n";
}

static inline void test_lfsr_hash_coverage_1() {
    io_u::io_utils io;
    lfsr_hash::gens g;
    constexpr int M = 64*4;
    std::array<uint8_t, M> buff{};
    std::cout << "Wait for LFSR 32-bit hashes coverage test 1.1...\n";
    std::set<lfsr8::u32> hashes;
    for (int i=0; i<256; i++) {
        const uint8_t x = i;
        buff[0] = x;
        hashes.insert(lfsr_hash::hash32<M>(g, buff.data(), {1 % M, 1, 1}));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / 256.) << "%" << std::endl;
    //
    for (int i=0; i<256*256; i++) {
        const lfsr8::u16 x = i;
        io.copy_to_mem_16(x, buff.data(), buff.size()); // LE guaranteed
        hashes.insert(lfsr_hash::hash32<M>(g, buff.data(), {2 % M, 2, 2}));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256.)) << "%" << std::endl;
    //
    for (int i=0; i<256*256*256; i++) {
        const lfsr8::u32 x = i;
        buff[0] = x % 256u;
        buff[1] = (x >> 8) % 256u;
        buff[2] = (x >> 16) % 256u;
        hashes.insert(lfsr_hash::hash32<M>(g, buff.data(), {3 % M, 3, 3}));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256. + 65536.*256.)) << "%" << std::endl;
    std::cout << " All Ok! Completed.\n";
}

static inline void test_lfsr_hash_coverage_2() {
    io_u::io_utils io;
    lfsr_hash::gens g;
    constexpr int M = 64*4;
    std::array<uint8_t, M> buff{};
    std::cout << "Wait for LFSR 64-bit hashes coverage test 1.2...\n";
    std::set<lfsr8::u64> hashes;
    for (int i=0; i<256; i++) {
        const uint8_t x = i;
        buff[0] = x;
        hashes.insert(lfsr_hash::hash64<M>(g, buff.data(), {1 % M, 1, 1}));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / 256.) << "%" << std::endl;
    //
    for (int i=0; i<256*256; i++) {
        const lfsr8::u16 x = i;
        io.copy_to_mem_16(x, buff.data(), buff.size()); // LE guaranteed
        hashes.insert(lfsr_hash::hash64<M>(g, buff.data(), {2 % M, 2, 2}));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256.)) << "%" << std::endl;
    //
    for (int i=0; i<256*256*256; i++) {
        const lfsr8::u32 x = i;
        buff[0] = x % 256u;
        buff[1] = (x >> 8) % 256u;
        buff[2] = (x >> 16) % 256u;
        hashes.insert(lfsr_hash::hash64<M>(g, buff.data(), {3 % M, 3, 3}));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256. + 65536.*256.)) << "%" << std::endl;
    std::cout << " All Ok! Completed.\n";
}
