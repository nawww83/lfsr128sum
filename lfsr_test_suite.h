#pragma once

#include <iostream>
#include <vector>
#include <string_view>
#include <algorithm>
#include <format> // C++20

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

    // Параметры Golden Hash
    // Вектор намеренно включает слитные блоки, 16-байтовый остаток и
    // побайтовый хвост, чтобы эталон фиксировал всю цепочку обработки.
    static constexpr size_t GOLDEN_M = 200;
    static constexpr lfsr_hash::u128 GOLDEN_EXPECTED
    {
        9682266398336795185ull, 1291011280817091015ull
    };

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

        test_fused_blocks();
        test_simd_step_equivalency();

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

    /**
     * @brief Скалярный эталон обработки входа (семантика next()).
     */
    static void process_input_scalar(lfsr_hash::gens &g, std::span<const std::byte> input)
    {
        const size_t n = input.size();
        if (n == 0) return;
        size_t i = 0;
        size_t j = n;

        for (; i + 16 <= n; i += 16)
        {
            for (int s = 0; s < 4; ++s)
            {
                const uint16_t lo = static_cast<uint16_t>(
                    static_cast<uint8_t>(input[i + 2 * s]) |
                    (static_cast<uint16_t>(static_cast<uint8_t>(input[i + 2 * s + 1])) << 8));
                const uint16_t hi = static_cast<uint16_t>(
                    static_cast<uint8_t>(input[i + 8 + 2 * s]) |
                    (static_cast<uint16_t>(static_cast<uint8_t>(input[i + 8 + 2 * s + 1])) << 8));

                g.g_251x4.next(lo, hi);
                g.g_241x4.next(static_cast<uint16_t>(~lo), static_cast<uint16_t>(~hi));
            }
        }

        while (i < j)
        {
            const uint16_t val_l = static_cast<uint8_t>(input[i]);
            const uint16_t val_r = (i + 1 < j) ? static_cast<uint8_t>(input[j - 1]) : 0;

            g.g_251x4.next(val_l, val_r);
            g.g_241x4.next(static_cast<uint16_t>(~val_r), static_cast<uint16_t>(~val_l));

            i++;
            if (j > i)
                j--;
        }

        const uint8_t first = static_cast<uint8_t>(input[0]);
        const uint16_t x = static_cast<uint16_t>(first | (first << 8));
        g.g_251x4.next(x, x);
        g.g_241x4.next(x, x);
    }

    /**
     * @brief Проверка слитного AVX2-прохода против скалярного эталона.
     * @details Покрывает длины, кратные и некратные размеру блока, а также
     * все ветки хвоста (16-байтовый остаток и побайтовый остаток).
     */
    void test_fused_blocks()
    {
        static const size_t lens[] = {
            1, 2, 3, 15, 16, 17, 31, 47, 63, 64, 65, 79, 127, 128, 129,
            191, 192, 255, 256, 1000, 4096, 65535, 65536, 65537};

        bool ok = true;
        std::string fail_len = "-";

        for (size_t len : lens)
        {
            std::vector<uint8_t> data(len);
            uint32_t x = 123456789u + static_cast<uint32_t>(len);
            for (size_t k = 0; k < len; ++k)
            {
                x = x * 1664525u + 1013904223u;
                data[k] = static_cast<uint8_t>(x >> 24);
            }

            lfsr_hash::gens g_fast, g_ref;
            g_fast.reset();
            g_ref.reset();
            g_fast.add_salt({5, 7, 9});
            g_ref.add_salt({5, 7, 9});

            const auto view = std::as_bytes(std::span<const uint8_t>(data));
            g_fast.process_input(view);
            process_input_scalar(g_ref, view);

            if (!(g_fast.g_251x4.get_state() == g_ref.g_251x4.get_state()) ||
                !(g_fast.g_241x4.get_state() == g_ref.g_241x4.get_state()))
            {
                ok = false;
                fail_len = std::to_string(len);
                break;
            }
        }

        report("Fused AVX2 vs scalar reference", ok,
               "all lengths match",
               ok ? "all lengths match" : ("mismatch at len " + fail_len));
    }

    /**
     * @brief Проверка одиночного SIMD-шага против скалярного next().
     * @details Важно, что состояние уже не единичное: иначе обратная связь
     * равна нулю и переполнение входа не проявляется.
     */
    void test_simd_step_equivalency()
    {
        lfsr_hash::gens g_simd, g_scalar;
        g_simd.reset();
        g_scalar.reset();
        for (int i = 0; i < 37; ++i)
        {
            g_simd.g_251x4.next(i, i + 3);
            g_scalar.g_251x4.next(i, i + 3);
        }

        uint32_t x = 1;
        bool ok = true;
        for (int i = 0; i < 4096; ++i)
        {
            x = x * 1664525u + 1013904223u;
            const uint16_t inp1 = static_cast<uint16_t>(x);
            const uint16_t inp2 = static_cast<uint16_t>(x >> 16);

            g_simd.g_251x4.next_simd(inp1, inp2);
            g_scalar.g_251x4.next(inp1, inp2);

            if (!(g_simd.g_251x4.get_state() == g_scalar.g_251x4.get_state()))
            {
                ok = false;
                break;
            }
        }

        report("SIMD step vs scalar next", ok, "match", ok ? "match" : "mismatch");
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

        const size_t total_attempts = 256 + 65536;

        std::vector<HashT> hashes;
        hashes.reserve(total_attempts);

        auto fill_hashes = [&](int count, int salt_val)
        {
            auto data_view = std::as_bytes(std::span(buff));
            g.reset();
            g.add_salt({salt_val % M, (uint8_t)salt_val, (uint8_t)salt_val});
            const auto& ref_state_241 = g.g_241x4.get_state();
            const auto& ref_state_251 = g.g_251x4.get_state();

            for (int i = 0; i < count; ++i)
            {
                io.copy_to_mem(i, buff);
                g.g_241x4.set_state(ref_state_241);
                g.g_251x4.set_state(ref_state_251);
                if constexpr (sizeof(HashT) == 4)
                    hashes.push_back(lfsr_hash::hash64(g, data_view)); // Автопреобразование в 32 бита из 64 битов.
                else
                    hashes.push_back(lfsr_hash::hash64(g, data_view));
            }
        };

        fill_hashes(256, 1);
        fill_hashes(65536, 2);

        // Считаем уникальные элементы через сортировку
        std::sort(hashes.begin(), hashes.end());
        auto last = std::unique(hashes.begin(), hashes.end());
        size_t unique_count = std::distance(hashes.begin(), last);

        // Форматируем результат через std::format (C++20)
        double coverage_pct = (static_cast<double>(unique_count) * 100.0) / total_attempts;

        // Для 32-битного хэша на 65792 входах по парадоксу дней рождения
        // ожидается ~0.5 коллизии, поэтому допускаем небольшой запас.
        // Для 64-битного хэша коллизии практически невозможны.
        const size_t allowed_collisions = (sizeof(HashT) == 4) ? 8 : 0;
        bool ok = (total_attempts - unique_count <= allowed_collisions);

        // Точный вывод: "100.00% (65792)"
        std::string actual_str = std::format("{:.2f}% ({} unique)", coverage_pct, unique_count);
        std::string expected_str = (allowed_collisions > 0)
            ? std::format(">= {} unique (допуск {})", total_attempts - allowed_collisions, allowed_collisions)
            : std::to_string(total_attempts);

        report(label, ok, expected_str, actual_str);
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
            buff[i] = static_cast<uint8_t>(i * 97 + 13);

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
