#pragma once

#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include "lfsr_hash.h"

#include <chrono>

namespace timer_n
{

    class Timer
    {
        std::chrono::time_point<std::chrono::steady_clock> mBegin;

    public:
        Timer() : mBegin(std::chrono::steady_clock::now()) {};
        double elapsed_ns()
        {
            auto mEnd = std::chrono::steady_clock::now();
            return std::chrono::duration_cast<std::chrono::nanoseconds>(mEnd - mBegin).count();
        }
        void reset()
        {
            mBegin = std::chrono::steady_clock::now();
        }
    };

}

template <size_t N>
double bench(timer_n::Timer &timer, lfsr_hash::gens &g, const std::vector<uint8_t> &data)
{
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
    if (result.first == 0xDEADBEEF)
        std::cout << " "; // Маловероятное условие, мешающее удалению кода
    return perf;
}

inline void lfsr_hash_benchmark()
{
    std::cout << "Wait for 128-bit LFSR hash benchmark...\n";
    std::cout << std::flush;
    timer_n::Timer timer;
    lfsr_hash::gens g;
    constexpr size_t N = 4 * 1024 * 1024;
    auto v = std::vector<uint8_t>(N);
    std::cout << "Input array of " << N << " bytes is allocated.\n";
    std::vector<double> perfs{};
    {
        // WarmUp
        const auto perf = bench<N>(timer, g, v);
    }
    constexpr size_t M = 128;
    for (int i = 0; i < M; ++i)
    {
        const auto perf = bench<N>(timer, g, v);
        if (std::isfinite(perf))
            perfs.push_back(perf);
    }
    std::sort(perfs.begin(), perfs.end());
    std::cout << "128-bit LFSR hash median performance: " << perfs[perfs.size() / 2] << " MB/s\n";
    std::cout << " All Ok! Completed.\n";
}

inline void check_hash(const std::string &version)
{
    lfsr_hash::gens g;
    constexpr int M = 64;
    std::array<uint8_t, M> buff;
    for (int i = 0; i < M; i++)
    {
        buff[i] = i;
    }
    g.reset();
    g.add_salt({8, 8, 8});
    lfsr_hash::u128 hash = lfsr_hash::hash128(g, std::as_bytes(std::span(buff)));
    if (version == "v2.1-simd") {
        assert((hash == lfsr_hash::u128{4355867762154551223ull, 17377028108451251887ull})); // Фиксация алгоритма, v.2.1.
        std::cout << "Golden hash: " << hash.first << ", " << hash.second << ", Ok\n";
    }
}

inline void test_simd_consistency()
{
    lfsr_hash::gens generator;
    generator.reset(); // Установит исходное состояние {1,0,0,0, 1,0,0,0}

    // Сохраняем начальное состояние для сравнения
    auto s0 = generator.g_251x4.get_state();

    // Используем явные типы и приведение
    const uint16_t val_input_1 = static_cast<uint16_t>(123 % 251);
    const uint16_t val_input_2 = static_cast<uint16_t>(200 % 241);

    // Шаг вперед и назад
    generator.g_251x4.next_simd(val_input_1, val_input_2);
    generator.g_251x4.back_simd(val_input_1, val_input_2);

    auto s_final = generator.g_251x4.get_state();

    bool passed = true;
    int idx = -1;
    for (int i = 0; i < 8; ++i)
    {
        if (s0[i] != s_final[i])
        {
            passed = false;
            idx = i;
            break;
        }
    }

    if (passed)
    {
        std::cout << "Consistency Test (next_simd <-> simd_back): PASSED" << std::endl;
    }
    else
    {
        std::cout << "Consistency Test: FAILED (Math mismatch at index " << idx << ")" << std::endl;
    }
}

inline void test_block_simd_consistency()
{
    lfsr_hash::gens generator;
    generator.reset(); // Устанавливаем начальное состояние {1,0,0,0, 1,0,0,0}

    // 1. Сохраняем исходное состояние (State_0)
    auto s_start = generator.g_251x4.get_state();

    // 2. Генерируем тестовый 128-битный блок данных (Data)
    // Используем разные значения для каждого байта, чтобы проверить диффузию
    alignas(16) uint8_t raw_data[16];
    for (int i = 0; i < 16; ++i)
    {
        raw_data[i] = static_cast<uint8_t>((i * 23 + 7) % 251);
    }
    __m128i data_block = _mm_loadu_si128(reinterpret_cast<const __m128i *>(raw_data));

    // 3. Выполняем шаг вперед на 16 байт сразу
    generator.g_251x4.next_simd_block(data_block);

    // 4. Выполняем шаг назад на те же 16 байт
    generator.g_251x4.back_simd_block(data_block);

    // 5. Проверка результата (State_0 == State_final)
    auto s_final = generator.g_251x4.get_state();
    bool passed = true;
    for (int i = 0; i < 8; ++i)
    {
        if (s_start[i] != s_final[i])
        {
            passed = false;
            break;
        }
    }

    std::cout << "--- Block SIMD Consistency Test (128-bit) ---" << std::endl;
    if (passed)
    {
        std::cout << "Status: PASSED" << std::endl;
    }
    else
    {
        std::cout << "Status: FAILED!" << std::endl;
        std::cout << "Expected: ";
        for (int i = 0; i < 8; ++i)
            std::cout << (int)s_start[i] << " ";
        std::cout << "\nActual:   ";
        for (int i = 0; i < 8; ++i)
            std::cout << (int)s_final[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "---------------------------------------------" << std::endl;
}

inline void test_lfsr_hash_coverage_1()
{
    io_u::io_utils io;
    lfsr_hash::gens g;
    constexpr int M = 64 * 4;
    std::array<uint8_t, M> buff;
    buff.fill(0);
    std::cout << "Wait for LFSR 32-bit hashes coverage test 1.1...\n";
    std::unordered_set<lfsr8::u32> hashes;
    for (int i = 0; i < 256; i++)
    {
        const uint8_t x = i;
        io.copy_to_mem(x, buff);
        g.reset();
        g.add_salt({1 % M, 1, 1});
        hashes.insert(lfsr_hash::hash32(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / 256.) << "%" << std::endl;
    for (int i = 0; i < 256 * 256; i++)
    {
        const lfsr8::u16 x = i;
        io.copy_to_mem(x, buff);
        g.reset();
        g.add_salt({2 % M, 2, 2});
        hashes.insert(lfsr_hash::hash32(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256.)) << "%" << std::endl;
    std::cout << " All Ok! Completed.\n";
}

inline void test_lfsr_hash_coverage_2()
{
    io_u::io_utils io;
    lfsr_hash::gens g;
    constexpr int M = 64 * 4;
    std::array<uint8_t, M> buff;
    buff.fill(0);
    std::cout << "Wait for LFSR 64-bit hashes coverage test 1.2...\n";
    std::unordered_set<lfsr8::u64> hashes;
    for (int i = 0; i < 256; i++)
    {
        const uint8_t x = i;
        io.copy_to_mem(x, buff);
        g.reset();
        g.add_salt({1 % M, 1, 1});
        hashes.insert(lfsr_hash::hash64(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / 256.) << "%" << std::endl;
    for (int i = 0; i < 256 * 256; i++)
    {
        const lfsr8::u16 x = i;
        io.copy_to_mem(x, buff);
        g.reset();
        g.add_salt({2 % M, 2, 2});
        hashes.insert(lfsr_hash::hash64(g, std::as_bytes(std::span(buff))));
    }
    std::cout << hashes.size() << ", coverage: " << (hashes.size() * 100. / (65536. + 256.)) << "%" << std::endl;
    std::cout << " All Ok! Completed.\n";
}
