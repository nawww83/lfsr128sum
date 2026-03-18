#pragma once

#include <iostream>
#include <vector>
#include <string_view>
#include <unordered_set>
#include <iomanip>

#include "timer.h"
#include "hardcores.h"
#include "lfsr_hash.h"
#include "io_utils.h"
#include "version.h"

// В CMakeLists.txt мы получили эти макросы из шаблона .in
#define VERSION_INT (PROJECT_VERSION_MAJOR_INT * 10000 + PROJECT_VERSION_MINOR_INT * 100 + PROJECT_VERSION_PATCH_INT)

class LFSRTestSuite
{
    int total = 0;
    int passed = 0;

    // ANSI-коды для красоты в терминале
    const std::string_view RESET = "\033[0m";
    const std::string_view RED = "\033[31m";
    const std::string_view GREEN = "\033[32m";
    const std::string_view CYAN = "\033[36m";
    const std::string_view GOLD = "\033[33m";

    // Параметры Golden Hash (v.2.1.0)
    static constexpr size_t GOLDEN_M = 64;
    static constexpr lfsr_hash::u128 GOLDEN_EXPECTED{
        17892976477579333464ull, 4582246380472290850ull};

    void report(std::string_view name, bool ok, std::string_view expected = "", std::string_view actual = "")
    {
        total++;
        if (ok)
            passed++;

        std::cout << "[" << (ok ? GREEN : RED) << (ok ? "OK  " : "FAIL") << RESET << "] " << name;
        if (!actual.empty())
            std::cout << " -> " << GOLD << actual << RESET;
        std::cout << "\n";

        if (!ok && !expected.empty())
        {
            std::cout << "      " << CYAN << "Expected: " << RESET << expected << "\n";
            std::cout << "      " << RED << "Actual:   " << RESET << actual << "\n";
        }
    }

    // Превращает состояние {1, 0, 2...} в строку "1 0 2..."
    std::string state_to_string(const auto &state)
    {
        std::string res;
        for (auto v : state)
        {
            res += std::to_string(static_cast<int>(v)) + " ";
        }
        return res;
    }

public:
    void run_all()
    {
        std::cout << CYAN << "=== Starting LFSR Test Suite ===" << RESET << "\n\n";

        test_scalar_and_simd_equivalency();
        test_scalar_and_simd_equivalency_overflowed_input();

        test_simd_consistency();
        test_block_simd_consistency();
        test_simd_consistency_overflowed_input();
        test_simd_block_consistency_overflowed_input();

        run_coverage_test<lfsr8::u32>("32-bit Hash Coverage");
        run_coverage_test<lfsr8::u64>("64-bit Hash Coverage");

        run_period_check();

        test_golden_hash();

        bool all_ok = (passed == total);
        std::cout << "\n"
                  << (all_ok ? GREEN : RED) << "=== SUMMARY: "
                  << passed << "/" << total << " tests passed ===" << RESET << "\n";
    }

    void test_scalar_and_simd_equivalency()
    {
        lfsr_hash::gens g_simd, g_scalar;
        g_simd.reset();
        g_scalar.reset();

        alignas(16) uint8_t raw[16];
        for (int i = 0; i < 16; ++i)
            raw[i] = (i * 23 + 7) % 251;
        __m128i block = _mm_load_si128((__m128i *)raw);

        g_simd.g_251x4.next_simd_block(block);

        alignas(16) uint16_t vals[8];
        _mm_store_si128((__m128i *)vals, block);
        for (int i = 0; i < 4; ++i)
            g_scalar.g_251x4.next(vals[i], vals[i + 4]);

        auto s_simd = g_simd.g_251x4.get_state();
        auto s_scalar = g_scalar.g_251x4.get_state();
        bool ok = (s_simd == s_scalar);

        report("Scalar vs SIMD Equivalency", ok,
               state_to_string(s_scalar),
               state_to_string(s_simd));
    }

    void test_scalar_and_simd_equivalency_overflowed_input()
    {
        lfsr_hash::gens g_simd, g_scalar;
        g_simd.reset();
        g_scalar.reset();

        alignas(16) uint8_t raw[16];
        for (int i = 0; i < 16; ++i)
            raw[i] = 255 - i;
        __m128i block = _mm_load_si128((__m128i *)raw);

        g_simd.g_251x4.next_simd_block(block);

        alignas(16) uint16_t vals[8];
        _mm_store_si128((__m128i *)vals, block);
        for (int i = 0; i < 4; ++i)
            g_scalar.g_251x4.next(vals[i], vals[i + 4]);

        auto s_simd = g_simd.g_251x4.get_state();
        auto s_scalar = g_scalar.g_251x4.get_state();
        bool ok = (s_simd == s_scalar);

        report("Scalar vs SIMD Equivalency for overflowed input", ok,
               state_to_string(s_scalar),
               state_to_string(s_simd));
    }

    void test_simd_consistency_overflowed_input()
    {
        lfsr_hash::gens g_simd, g_scalar;
        g_simd.reset();
        g_scalar.reset();

        alignas(16) uint8_t raw[16];
        for (int i = 0; i < 16; ++i)
            raw[i] = 255 - i;
        __m128i block = _mm_load_si128((__m128i *)raw);

        alignas(16) uint16_t vals[8];
        _mm_store_si128((__m128i *)vals, block);
        for (int i = 0; i < 4; ++i)
            g_scalar.g_251x4.next(vals[i], vals[i + 4]);
        for (int i = 4 - 1; i >= 0; --i)
            g_scalar.g_251x4.back(vals[i], vals[i + 4]);

        for (int i = 0; i < 4; ++i)
            g_simd.g_251x4.next_simd(vals[i], vals[i + 4]);
        for (int i = 4 - 1; i >= 0; --i)
            g_simd.g_251x4.back_simd(vals[i], vals[i + 4]);

        auto s_simd = g_simd.g_251x4.get_state();
        auto s_scalar = g_scalar.g_251x4.get_state();
        bool ok = (s_simd == s_scalar) && (s_simd == lfsr8::u16x8{1, 0, 0, 0, 1, 0, 0, 0});

        report("SIMD Next/Back Consistency for overflowed input", ok,
               state_to_string(s_scalar),
               state_to_string(s_simd));
    }

    void test_simd_block_consistency_overflowed_input()
    {
        lfsr_hash::gens g_simd, g_scalar;
        g_simd.reset();
        g_scalar.reset();

        alignas(16) uint8_t raw[16];
        for (int i = 0; i < 16; ++i)
            raw[i] = 255 - i;
        __m128i block = _mm_load_si128((__m128i *)raw);

        g_simd.g_251x4.next_simd_block(block);
        g_simd.g_251x4.back_simd_block(block);

        alignas(16) uint16_t vals[8];
        _mm_store_si128((__m128i *)vals, block);
        for (int i = 0; i < 4; ++i)
            g_scalar.g_251x4.next(vals[i], vals[i + 4]);
        for (int i = 4 - 1; i >= 0; --i)
            g_scalar.g_251x4.back(vals[i], vals[i + 4]);

        auto s_simd = g_simd.g_251x4.get_state();
        auto s_scalar = g_scalar.g_251x4.get_state();
        bool ok = (s_simd == s_scalar) && (s_simd == lfsr8::u16x8{1, 0, 0, 0, 1, 0, 0, 0});

        report("Block SIMD Consistency for overflowed input", ok,
               state_to_string(s_scalar),
               state_to_string(s_simd));
    }

    void test_simd_consistency()
    {
        lfsr_hash::gens gen;
        gen.reset();
        const auto s0 = gen.g_251x4.get_state();

        for (int i = 0; i < 8; ++i)
            gen.g_251x4.next_simd(i, i + 1);
        for (int i = 8 - 1; i >= 0; --i)
            gen.g_251x4.back_simd(i, i + 1);

        const auto s_final = gen.g_251x4.get_state();

        report("SIMD Next/Back Consistency",
               s_final == s0,
               state_to_string(s0),
               state_to_string(s_final));
    }

    void test_block_simd_consistency()
    {
        lfsr_hash::gens gen;
        gen.reset();
        const auto s0 = gen.g_251x4.get_state();

        alignas(16) uint8_t raw[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
        __m128i block = _mm_load_si128(reinterpret_cast<const __m128i *>(raw));

        gen.g_251x4.next_simd_block(block);
        gen.g_251x4.next_simd_block(block);
        gen.g_251x4.next_simd_block(block);
        gen.g_251x4.next_simd_block(block);
        gen.g_251x4.next_simd_block(block);
        gen.g_251x4.next_simd_block(block);
        gen.g_251x4.next_simd_block(block);
        gen.g_251x4.next_simd_block(block);

        gen.g_251x4.back_simd_block(block);
        gen.g_251x4.back_simd_block(block);
        gen.g_251x4.back_simd_block(block);
        gen.g_251x4.back_simd_block(block);
        gen.g_251x4.back_simd_block(block);
        gen.g_251x4.back_simd_block(block);
        gen.g_251x4.back_simd_block(block);
        gen.g_251x4.back_simd_block(block);

        const auto s_final = gen.g_251x4.get_state();

        report("Block SIMD Consistency",
               s_final == s0,
               state_to_string(s0),
               state_to_string(s_final));
    }

    template <typename HashT>
    void run_coverage_test(std::string_view label)
    {
        io_u::io_utils io;
        lfsr_hash::gens g;
        constexpr int M = 64 * 4;
        std::array<uint8_t, M> buff;
        buff.fill(0);

        std::unordered_set<HashT> hashes;
        hashes.reserve(256 + 65536);

        // Лямбда просто наполняет множество
        auto fill_hashes = [&](int count, int salt_val)
        {
            for (int i = 0; i < count; ++i)
            {
                io.copy_to_mem(i, buff);
                g.reset();
                g.add_salt({salt_val % M, (uint8_t)salt_val, (uint8_t)salt_val});

                if constexpr (sizeof(HashT) == 4)
                    hashes.insert(lfsr_hash::hash32(g, std::as_bytes(std::span(buff))));
                else
                    hashes.insert(lfsr_hash::hash64(g, std::as_bytes(std::span(buff))));
            }
        };

        // Выполняем два прохода (накопительно)
        fill_hashes(256, 1);
        fill_hashes(65536, 2);

        // Итоговый расчет покрытия
        const size_t total_attempts = 256 + 65536;
        const double coverage_pct = (static_cast<double>(hashes.size()) * 100.0) / total_attempts;

        bool ok = (hashes.size() == total_attempts);

        // Формируем строку результата: "100.00% (65792/65792)"
        std::string actual_str = std::to_string(coverage_pct).substr(0, 5) + "% (" + std::to_string(hashes.size()) + ")";

        report(label, ok, std::to_string(total_attempts), actual_str);
    }

    void run_period_check()
    {
        lfsr_hash::gens g;
        const auto &K1 = lfsr_hash::K1;
        const auto &K2 = lfsr_hash::K2;
        // Правый генератор имеет период T1 = p^(m-1) - 1 (малый период). У него есть "замороженные" состояния c*(K0, K0 + K1, K0 + K1 + K2, 1).
        // c - константа, [0, p-1].
        // Это так, потому что полином g(x) = (x - 1)v(x), где v(x) - полином (m-1)-степени с наибольшим периодом T1.
        // Полином (x-1) имеет длину цикла равную единице - он и порождает "замороженные" состояния.
        // У нас LFSR реализован в форме Галуа. Т.о., мы косвенно тестируем генератор на т.н. малый период.
        lfsr_hash::STATE lock_up_state1 = {0u, 0u, 0u, 0u, K2[4], (K2[4] + K2[5]) % 241u, (K2[4] + K2[5] + K2[6]) % 241u, 1u};
        g.g_241x4.set_state(lock_up_state1);
        g.g_241x4.next();
        bool ok1 = g.g_241x4.is_state_high(lock_up_state1);

        lfsr_hash::STATE lock_up_state2 = {0u, 0u, 0u, 0u, K1[4], (K1[4] + K1[5]) % 251u, (K1[4] + K1[5] + K1[6]) % 251u, 1u};
        g.g_251x4.set_state(lock_up_state2);
        g.g_251x4.next();
        bool ok2 = g.g_251x4.is_state_high(lock_up_state2);
        report("Period p^(m-1) test LFSR p = 241", ok1, state_to_string(lock_up_state1), state_to_string(g.g_241x4.get_state()));
        report("Period p^(m-1) test LFSR p = 251", ok2, state_to_string(lock_up_state2), state_to_string(g.g_251x4.get_state()));
    }

    void run_lfsr_benchmark()
    {
        hardcores::bind_to_core(1); // Стабилизация показаний.
        constexpr size_t N = 8 * 1024 * 1024;
        lfsr_hash::gens g;
        static std::vector<uint8_t> v(N, 0xAA);

        // Прогрев кэша.
        lfsr_hash::hash128(g, std::as_bytes(std::span(v)));

        timer_n::Timer timer;

        g.reset();
        g.add_salt({1, 2, 3});
        auto result = lfsr_hash::hash128(g, std::as_bytes(std::span(v)));
        hardcores::doNotOptimizeAway(result);
        const auto dt_ns = timer.elapsed_ns();

        double mb_s = (static_cast<double>(N) / (1024.0 * 1024.0)) / (dt_ns / 1e9);
        report("128-bit Performance", true, "", std::to_string((int)mb_s) + " MB/s");
    }

    void test_golden_hash()
    {
        lfsr_hash::gens g;
        std::array<uint8_t, GOLDEN_M> buff;
        for (int i = 0; i < GOLDEN_M; i++)
            buff[i] = static_cast<uint8_t>(i);

        g.reset();
        g.add_salt({8, 8, 8});
        lfsr_hash::u128 hash = lfsr_hash::hash128(g, std::as_bytes(std::span(buff)));

        // Проверка версии через макрос из CMake
        if constexpr (VERSION_INT >= 20000)
        {
            bool ok = (hash == GOLDEN_EXPECTED);

            auto to_str = [](const lfsr_hash::u128 &h)
            {
                return std::to_string(h.first) + ", " + std::to_string(h.second);
            };

            report("Golden Hash", ok, to_str(GOLDEN_EXPECTED), to_str(hash));
        }
        else
        {
            report("Golden Hash", true, "", "Skipped");
        }
    }
};
