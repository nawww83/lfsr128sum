#pragma once

#include <cmath>
#include <iostream>
#include <set>
#include <unordered_set>
#include <vector>
#include <algorithm>
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
    if (data.size() > N) 
    {
        return (std::numeric_limits<double>::infinity)();
    }
    timer.reset();
    g.reset();
    g.add_salt({1, 2, 3});
    lfsr_hash::u128 result = lfsr_hash::hash128(g, std::as_bytes(std::span(data)));
    const auto dt = timer.elapsed_ns();
    const double perf = (1000. * N) / dt;
    if (result.first == 0xDEADBEEF) std::cout << " "; // Маловероятное условие, мешающее удалению кода
    return perf;
}

inline void lfsr_hash_benchmark() 
{
    std::cout << "Wait for 128-bit LFSR hash benchmark...\n";
    std::cout << std::flush;
    timer_n::Timer timer;
    lfsr_hash::gens g;
    constexpr size_t N = 4*1024*1024;
    auto v = std::vector<uint8_t>(N);
    std::cout << "Input array of " << N << " bytes is allocated.\n";
    std::vector<double> perfs {};
    { 
        // WarmUp
        const auto perf = bench<N>(timer, g, v);
    }
    constexpr size_t M = 128;
    for (int i = 0; i < M; ++i) {
        const auto perf = bench<N>(timer, g, v);
        if (std::isfinite(perf))
            perfs.push_back( perf );
    }
    std::sort(perfs.begin(), perfs.end());
    std::cout << "128-bit LFSR hash median performance: " << perfs[perfs.size()/2] << " MB/s\n";
    std::cout << " All Ok! Completed.\n";
}

inline void test_lfsr_hash_coverage_1() {
    io_u::io_utils io;
    lfsr_hash::gens g;
    constexpr int M = 64*4;
    std::array<uint8_t, M> buff{};
    std::cout << "Wait for LFSR 32-bit hashes coverage test 1.1...\n";
    std::unordered_set<lfsr8::u32> hashes;
    // hashes.reserve(256ull*256ull*256ull);
    for (int i=0; i<256; i++) {
        const uint8_t x = i;
        buff[0] = x;
        g.reset();
        g.add_salt({1 % M, 1, 1});
        hashes.insert(lfsr_hash::hash32(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / 256.) << "%" << std::endl;
    for (int i=0; i<256*256; i++) {
        const lfsr8::u16 x = i;
        io.copy_to_mem(x, buff);
        g.reset();
        g.add_salt({2 % M, 2, 2});
        hashes.insert(lfsr_hash::hash32(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256.)) << "%" << std::endl;
    // std::cout << "You should wait a bit longer...\n" << std::flush;
    // for (int i=0; i<256*256*256; i++) {
    //     const lfsr8::u32 x = i;
    //     buff[0] = x % 256u;
    //     buff[1] = (x >> 8) % 256u;
    //     buff[2] = (x >> 16) % 256u;
    //     g.reset();
    //     g.add_salt({3 % M, 3, 3});
    //     hashes.insert(lfsr_hash::hash32(g, std::as_bytes(std::span(buff))));
    // }
    // std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256. + 65536.*256.)) << "%" << std::endl;
    std::cout << " All Ok! Completed.\n";
}

inline void test_lfsr_hash_coverage_2() {
    io_u::io_utils io;
    lfsr_hash::gens g;
    constexpr int M = 64*4;
    std::array<uint8_t, M> buff{};
    std::cout << "Wait for LFSR 64-bit hashes coverage test 1.2...\n";
    std::unordered_set<lfsr8::u64> hashes;
    // hashes.reserve(256ull*256ull*256ull);
    for (int i=0; i<256; i++) {
        const uint8_t x = i;
        buff[0] = x;
        g.reset();
        g.add_salt({1 % M, 1, 1});
        hashes.insert(lfsr_hash::hash64(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / 256.) << "%" << std::endl;
    for (int i=0; i<256*256; i++) {
        const lfsr8::u16 x = i;
        io.copy_to_mem(x, buff);
        g.reset();
        g.add_salt({2 % M, 2, 2});
        hashes.insert(lfsr_hash::hash64(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256.)) << "%" << std::endl;
    // std::cout << "You should wait a bit longer...\n" << std::flush;
    // for (int i=0; i<256*256*256; i++) {
    //     const lfsr8::u32 x = i;
    //     buff[0] = x % 256u;
    //     buff[1] = (x >> 8) % 256u;
    //     buff[2] = (x >> 16) % 256u;
    //     g.reset();
    //     g.add_salt({3 % M, 3, 3});
    //     hashes.insert(lfsr_hash::hash64(g, std::as_bytes(std::span(buff))));
    // }
    // std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256. + 65536.*256.)) << "%" << std::endl;
    std::cout << " All Ok! Completed.\n";
}
