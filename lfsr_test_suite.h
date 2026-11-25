#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <span>
#include <string>
#include <string_view>
#include <algorithm>
#include <cstdio>

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
    // Вектор намеренно включает два 512-байтных слитных блока (то есть
    // два перекрёстных перемешивания полос фолда), 16-байтовый остаток и
    // побайтовый хвост, чтобы эталон фиксировал всю цепочку обработки.
    static constexpr size_t GOLDEN_M = 1037;
    static constexpr lfsr_hash::u128 GOLDEN_EXPECTED
    {
        15264202950503091111ull, 9712154424917034579ull
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

    /// Форматирование без <format> (GCC 11 и другие без std::format).
    static std::string fixed2(double v)
    {
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.2f", v);
        return std::string(buf);
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
        test_fused_salt_matches_power_by();
        test_simd_step_equivalency();

        run_coverage_test<lfsr8::u32>("32-bit Hash Coverage");
        run_coverage_test<lfsr8::u64>("64-bit Hash Coverage");

        run_period_check();

        test_golden_hash();
        test_start_prototype_matches_salt();
        test_salt_injectivity();
        test_fold_lane_mixing();

        test_mac_known_answer();
        test_mac_key_avalanche();
        test_mac_streaming_matches_oneshot();
        test_mac_key_lengths();
        test_mac_key_length_binding();
        test_linear_injection_attack();
        test_streaming_consistency();
        test_length_binding();
        test_zero_salt_does_not_reset();

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
     * @brief Линейная загрузка входа так же, как в hasher: слитные блоки
     * по 512 байт + остаток через process_remainder().
     */
    static void ingest_fast(lfsr_hash::gens &g, std::span<const std::byte> input)
    {
        constexpr size_t B = lfsr8::LFSR_fused_2x4<251, 241>::BLOCK_BYTES;
        const size_t nblocks = input.size() / B;
        if (nblocks > 0)
        {
            lfsr_hash::gens::fused_engine().process_blocks(
                input.data(), nblocks,
                g.g_251x4.state_ref(), g.g_241x4.state_ref());
        }
        lfsr_hash::process_remainder(g, input.subspan(nblocks * B));
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
            ingest_fast(g_fast, view);
            process_input_scalar(g_ref, view);

            if (!(g_fast.g_251x4.get_state() == g_ref.g_251x4.get_state()) ||
                !(g_fast.g_241x4.get_state() == g_ref.g_241x4.get_state()))
            {
                ok = false;
                fail_len = std::to_string(len);
                break;
            }
        }

        report("Fused engine (AVX2/SSE) vs scalar reference", ok,
               "all lengths match",
               ok ? "all lengths match" : ("mismatch at len " + fail_len));
    }

    /**
     * @brief Слитное возведение в степень в `add_salt` обязано побитно
     * совпадать с раздельными `power_by` для обоих генераторов.
     * @details v4 применяет соль к состояниям 251-й и 241-й пар одним
     * 256-битным проходом (LFSR_fused_exp). Здесь та же соль применяется
     * по отдельности, и состояния сравниваются после каждой соли.
     */
    void test_fused_salt_matches_power_by()
    {
        lfsr_hash::gens fused, plain;
        fused.reset();
        plain.reset();

        std::array<uint8_t, 128> data;
        for (size_t i = 0; i < data.size(); ++i)
            data[i] = static_cast<uint8_t>(i * 71 + 13);
        const auto view = std::as_bytes(std::span(data));
        ingest_fast(fused, view);
        ingest_fast(plain, view);

        bool ok = true;
        const lfsr_hash::salt salts[] = {lfsr_hash::S0, lfsr_hash::S1,
                                         lfsr_hash::S2};
        for (const auto &s : salts)
        {
            fused.add_salt(s);

            plain.g_251x4.next_simd(s.s0);
            plain.g_241x4.next_simd(s.s1);
            plain.g_251x4.power_by(s.q);
            plain.g_241x4.power_by(s.q);

            if (!(fused.g_251x4.get_state() == plain.g_251x4.get_state()) ||
                !(fused.g_241x4.get_state() == plain.g_241x4.get_state()))
                ok = false;
        }

        report("Fused salt vs separate power_by", ok,
               "all salts match", ok ? "all salts match" : "MISMATCH");
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
        constexpr int M = 64 * 4;
        std::array<uint8_t, M> buff;
        buff.fill(0);

        const size_t total_attempts = 256 + 65536;

        std::vector<HashT> hashes;
        hashes.reserve(total_attempts);

        auto fill_hashes = [&](int count, uint64_t base)
        {
            for (int i = 0; i < count; ++i)
            {
                io.copy_to_mem(base + static_cast<uint64_t>(i), buff);
                const auto h = lfsr_hash::hash128(std::as_bytes(std::span(buff)));
                hashes.push_back(static_cast<HashT>(h.first));
            }
        };

        fill_hashes(256, 0);
        fill_hashes(65536, 256);

        // Считаем уникальные элементы через сортировку
        std::sort(hashes.begin(), hashes.end());
        auto last = std::unique(hashes.begin(), hashes.end());
        size_t unique_count = std::distance(hashes.begin(), last);

        double coverage_pct = (static_cast<double>(unique_count) * 100.0) / total_attempts;

        // Для 32-битного хэша на 65792 входах по парадоксу дней рождения
        // ожидается ~0.5 коллизии, поэтому допускаем небольшой запас.
        // Для 64-битного хэша коллизии практически невозможны.
        const size_t allowed_collisions = (sizeof(HashT) == 4) ? 8 : 0;
        bool ok = (total_attempts - unique_count <= allowed_collisions);

        // Точный вывод: "100.00% (65792)"
        std::string actual_str = fixed2(coverage_pct) + "% (" +
                                 std::to_string(unique_count) + " unique)";
        std::string expected_str = (allowed_collisions > 0)
            ? ">= " + std::to_string(total_attempts - allowed_collisions) +
                  " unique (допуск " + std::to_string(allowed_collisions) + ")"
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
        static std::vector<uint8_t> v(N, 0xAA);

        // Прогрев кэша.
        lfsr_hash::hash128(std::as_bytes(std::span(v)));

        timer_n::Timer timer;

        auto result = lfsr_hash::hash128(std::as_bytes(std::span(v)));
        hardcores::doNotOptimizeAway(result);
        const auto dt_ns = timer.elapsed_ns();

        double mb_s = (static_cast<double>(N) / (1024.0 * 1024.0)) / (dt_ns / 1e9);
        report("128-bit Performance", true, "", std::to_string((int)mb_s) + " MB/s");
    }

    void test_golden_hash()
    {
        std::array<uint8_t, GOLDEN_M> buff;
        for (int i = 0; i < GOLDEN_M; i++)
            buff[i] = static_cast<uint8_t>(i * 97 + 13);

        lfsr_hash::u128 hash = lfsr_hash::hash128(std::as_bytes(std::span(buff)));

        bool ok = (hash == GOLDEN_EXPECTED);
        auto to_str = [](const lfsr_hash::u128 &h)
        {
            return std::to_string(h.first) + ", " + std::to_string(h.second);
        };
        report("Golden Hash v4", ok, to_str(GOLDEN_EXPECTED), to_str(hash));
    }

    /**
     * @brief Предвычисленное стартовое состояние обязано совпадать
     * с единичным состоянием после соли S0.
     * @details `hasher::reset()` берёт готовый прототип вместо возведения
     * в степень; тест ловит расхождение, если константы или путь соли
     * поменяются, а прототип забудут пересчитать.
     */
    void test_start_prototype_matches_salt()
    {
        lfsr_hash::gens ref;
        ref.reset();
        ref.add_salt(lfsr_hash::S0);

        const auto &start = lfsr_hash::gens::start();
        const bool ok = (start.g_251x4.get_state() == ref.g_251x4.get_state()) &&
                        (start.g_241x4.get_state() == ref.g_241x4.get_state());
        report("Start prototype == unit + S0", ok, "equal", ok ? "equal" : "MISMATCH");
    }

    static std::string u128_str(const lfsr_hash::u128 &h)
    {
        return std::to_string(h.first) + ", " + std::to_string(h.second);
    }

    static std::array<uint8_t, 32> key_bytes_0_31()
    {
        std::array<uint8_t, 32> k{};
        for (int i = 0; i < 32; ++i)
            k[i] = static_cast<uint8_t>(i);
        return k;
    }

    /**
     * @brief KAT для keyed-режима (MAC).
     * @details Ключ + сообщение -> фиксированный тег. Любая правка схемы
     * (стартовое состояние, сид фолда, слова финализации, доменная метка)
     * обязана ломать эти векторы, иначе формат MAC меняется молча.
     */
    void test_mac_known_answer()
    {
        const auto k32 = key_bytes_0_31();
        const auto kc = lfsr_hash::make_key(std::as_bytes(std::span(k32)));
        const char abc[] = "abc";

        const lfsr_hash::u128 e1{4482865992510714747ull, 14702203991202938416ull};
        const auto h1 = lfsr_hash::mac128(kc, std::as_bytes(std::span(abc, 3)));
        report("MAC KAT: k32 + \"abc\"", h1 == e1, u128_str(e1), u128_str(h1));

        std::array<uint8_t, GOLDEN_M> g{};
        for (size_t i = 0; i < g.size(); ++i)
            g[i] = static_cast<uint8_t>(i * 97 + 13);
        const lfsr_hash::u128 e2{17548153652833865696ull, 6324480247405512134ull};
        const auto h2 = lfsr_hash::mac128(kc, std::as_bytes(std::span(g)));
        report("MAC KAT: k32 + 1037 B", h2 == e2, u128_str(e2), u128_str(h2));

        uint8_t k1[1] = {0x42};
        const auto kc1 = lfsr_hash::make_key(std::as_bytes(std::span(k1)));
        const lfsr_hash::u128 e3{12238301807564121565ull, 5018636314847259397ull};
        const auto h3 = lfsr_hash::mac128(kc1, std::as_bytes(std::span<const uint8_t>()));
        report("MAC KAT: k1 + empty", h3 == e3, u128_str(e3), u128_str(h3));
    }

    /**
     * @brief Один бит ключа меняет тег целиком, и MAC не равен unkeyed-хэшу.
     */
    void test_mac_key_avalanche()
    {
        auto k32 = key_bytes_0_31();
        const auto kc = lfsr_hash::make_key(std::as_bytes(std::span(k32)));
        auto k32b = k32;
        k32b[31] ^= 1;
        const auto kcb = lfsr_hash::make_key(std::as_bytes(std::span(k32b)));

        std::vector<uint8_t> msg(64, 0xA5);
        const auto h1 = lfsr_hash::mac128(kc, std::as_bytes(std::span(msg)));
        const auto h2 = lfsr_hash::mac128(kcb, std::as_bytes(std::span(msg)));
        const int d = static_cast<int>(std::popcount(h1.first ^ h2.first) +
                                       std::popcount(h1.second ^ h2.second));

        // Для 128-битного выхода ожидание 64 бита; границы 40..88 — это
        // примерно ±4 сигмы, то есть ловят «почти совпало», не флакая шум.
        report("MAC: 1-bit key change", h1 != h2 && d >= 40 && d <= 88, "40..88 bits",
               std::to_string(d) + " bits");

        const auto hu = lfsr_hash::hash128(std::as_bytes(std::span(msg)));
        report("MAC != unkeyed hash", h1 != hu, "different",
               h1 != hu ? "different" : "SAME");
    }

    /**
     * @brief Keyed-хэшер потоковый: результат не зависит от нарезки,
     * а reset() возвращает то же ключевое стартовое состояние.
     */
    void test_mac_streaming_matches_oneshot()
    {
        const auto k32 = key_bytes_0_31();
        const auto kc = lfsr_hash::make_key(std::as_bytes(std::span(k32)));

        std::array<uint8_t, 5000> buf{};
        for (size_t i = 0; i < buf.size(); ++i)
            buf[i] = static_cast<uint8_t>(i * 37 + 11);

        const auto one = lfsr_hash::mac128(kc, std::as_bytes(std::span(buf)));

        lfsr_hash::hasher h(kc);
        const size_t steps[] = {1, 7, 13, 64, 511, 512, 513, 1000};
        size_t p = 0;
        for (size_t s : steps)
        {
            const size_t take = (std::min)(s, buf.size() - p);
            h.update(std::as_bytes(std::span(buf.data() + p, take)));
            p += take;
        }
        h.update(std::as_bytes(std::span(buf.data() + p, buf.size() - p)));
        const auto chunked = h.finish();

        report("MAC: streaming == one-shot", chunked == one, u128_str(one), u128_str(chunked));

        // Тот же объект после reset() обязан дать тот же тег.
        h.reset();
        h.update(std::as_bytes(std::span(buf)));
        const auto again = h.finish();
        report("MAC: reset() keeps the key", again == one, u128_str(one), u128_str(again));
    }

    /**
     * @brief Допустимые длины ключа и отказ на границах.
     */
    void test_mac_key_lengths()
    {
        bool ok = true;
        std::string detail;

        for (size_t n : {1u, 16u, 32u, 64u})
        {
            std::array<uint8_t, lfsr_hash::KEY_MAX> k{};
            for (size_t i = 0; i < n; ++i)
                k[i] = static_cast<uint8_t>(i * 7 + 3);
            try
            {
                const auto kc = lfsr_hash::make_key(std::as_bytes(std::span(k.data(), n)));
                ok &= !(kc.start.s251 == lfsr_hash::STATE{});
            }
            catch (...)
            {
                ok = false;
                detail += "len " + std::to_string(n) + " threw; ";
            }
        }

        bool threw_empty = false, threw_long = false;
        try
        {
            lfsr_hash::make_key(std::span<const std::byte>());
        }
        catch (const std::invalid_argument &)
        {
            threw_empty = true;
        }
        std::array<uint8_t, lfsr_hash::KEY_MAX + 1> too_long{};
        try
        {
            lfsr_hash::make_key(std::as_bytes(std::span(too_long)));
        }
        catch (const std::invalid_argument &)
        {
            threw_long = true;
        }

        ok = ok && threw_empty && threw_long;
        if (!ok)
            detail += std::string(threw_empty ? "" : "empty accepted; ") +
                      (threw_long ? "" : "65 B accepted; ");
        report("MAC: key length 1..64", ok, "1..64 ok, 0/65 rejected",
               ok ? "1..64 ok, 0/65 rejected" : detail);
    }

    /**
     * @brief Длина ключа входит в доменную метку: "a" и "a\0" — разные ключи.
     * @details Ключ дополняется нулями до 64 байт (как в HMAC), поэтому без
     * длины в метке эти два ключа были бы неразличимы.
     */
    void test_mac_key_length_binding()
    {
        const auto k32 = key_bytes_0_31();
        const auto kc = lfsr_hash::make_key(std::as_bytes(std::span(k32)));

        const uint8_t a1[2] = {'a', 0};
        const auto h1 = lfsr_hash::mac128(kc, std::as_bytes(std::span(a1, 1)));
        const auto h2 = lfsr_hash::mac128(kc, std::as_bytes(std::span(a1, 2)));

        report("MAC: \"a\" != \"a\\0\"", h1 != h2, "different",
               h1 != h2 ? "different" : "SAME");
    }

    /**
     * @brief Показатели солей обязаны быть взаимно просты со всеми периодами.
     * @details При gcd(q, T) > 1 отображение x -> x^q имеет ядро размера
     * gcd(q, T): образ сжимается, и вместе с ним падает эффективное
     * пространство состояний. Инвариант закреплён static_assert в заголовке,
     * здесь он просто виден в списке тестов.
     */
    void test_salt_injectivity()
    {
        const lfsr_hash::salt salts[] = {lfsr_hash::S0, lfsr_hash::S1,
                                         lfsr_hash::S2};
        bool ok = true;
        for (const auto &s : salts)
            ok &= lfsr_hash::salt_is_injective(s.q);
        report("Salt exponents coprime with periods", ok, "injective",
               ok ? "injective" : "NOT INJECTIVE");
    }

    /**
     * @brief Перекрёстная связь полос нелинейного фолда.
     * @details Различие в первом 32-байтном блоке (полоса 0) после полного
     * 512-байтного цикла (16 блоков, четыре прохода по полосам) обязано дойти
     * до всех четырёх полос: иначе след остаётся в одной цепочке.
     */
    void test_fold_lane_mixing()
    {
        std::array<uint8_t, 512> a;
        for (size_t i = 0; i < a.size(); ++i)
            a[i] = static_cast<uint8_t>(i * 31 + 7);
        auto b = a;
        b[0] ^= 1;

        lfsr_hash::fold_acc fa, fb;
        fa.reset();
        fb.reset();
        fa.update(a.data(), a.size());
        fb.update(b.data(), b.size());
        fa.flush();
        fb.flush();

        int diff = 0;
        for (int i = 0; i < 4; ++i)
        {
            lfsr8::u32 x[8], y[8];
            lfsr_hash::fold_store(fa.lane(i), x);
            lfsr_hash::fold_store(fb.lane(i), y);
            bool d = false;
            for (int j = 0; j < 8; ++j)
                d |= (x[j] != y[j]);
            diff += d ? 1 : 0;
        }
        report("Fold: cross-lane mixing", diff == 4, "4/4 lanes differ",
               std::to_string(diff) + "/4 lanes differ");
    }

    /**
     * @brief Регрессия на линейную инъекцию.
     * @details Без нелинейного слоя состояние зависело бы только от word mod p, поэтому правка
     * одного 16-битного слова на lcm(251,241) = 60491 не меняла хэш вообще
     * (бесплатная коллизия). В v4 нелинейный фолд обязан это поймать.
     */
    void test_linear_injection_attack()
    {
        constexpr uint32_t D = 251u * 241u; // 60491
        constexpr size_t N = 4096;

        auto make = [](std::array<uint8_t, N> &x, size_t off, uint32_t word)
        {
            for (size_t i = 0; i < N; ++i)
                x[i] = static_cast<uint8_t>(i * 37 + 11);
            x[off] = static_cast<uint8_t>(word & 0xFF);
            x[off + 1] = static_cast<uint8_t>(word >> 8);
        };

        bool ok = true;
        for (size_t off : {size_t(2), size_t(1000), size_t(2000), size_t(4094)})
        {
            std::array<uint8_t, N> a, b;
            make(a, off, 16);
            make(b, off, 16 + D);
            ok &= (lfsr_hash::hash128(std::as_bytes(std::span(a))) !=
                   lfsr_hash::hash128(std::as_bytes(std::span(b))));
        }
        report("v4: linear injection attack is detected", ok,
               "different hashes", ok ? "different hashes" : "COLLISION");
    }

    /// Потоковый режим обязан совпадать с одноразовым при любом размере кусков.
    void test_streaming_consistency()
    {
        constexpr size_t N = 100 * 1024 + 7;
        std::vector<uint8_t> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = static_cast<uint8_t>(i * 131 + 7);

        const auto one = lfsr_hash::hash128(std::as_bytes(std::span(data)));

        bool ok = true;
        for (size_t chunk : {size_t(1), size_t(7), size_t(32), size_t(511),
                             size_t(512), size_t(513), size_t(4096), size_t(65536)})
        {
            lfsr_hash::hasher h;
            for (size_t i = 0; i < N; i += chunk)
            {
                const size_t take = (std::min)(chunk, N - i);
                h.update(std::as_bytes(std::span(data.data() + i, take)));
            }
            if (h.finish() != one)
                ok = false;
        }
        report("v4: streaming == one-shot", ok, "equal", ok ? "equal" : "MISMATCH");
    }


    /**
     * @brief Регрессия: соль с q = 0 не должна стирать состояние.
     * @details Возведение в степень с нулём — тождество, а не x^0 = 1:
     * состояние не должно схлопываться в единичное.
     */
    void test_zero_salt_does_not_reset()
    {
        lfsr_hash::gens g;
        g.reset();

        // Уводим состояние от единичного.
        std::array<std::byte, 64> data{};
        ingest_fast(g, std::span<const std::byte>(data));

        const auto s251_before = g.g_251x4.get_state();
        const auto s241_before = g.g_241x4.get_state();

        g.g_251x4.power_by(0);
        g.g_241x4.power_by(0);
        const bool power_ok = (g.g_251x4.get_state() == s251_before) &&
                              (g.g_241x4.get_state() == s241_before);

        // Соль с q = 0: шаг с входом применяется, но состояние не обнуляется.
        g.add_salt({0, 2, 3});
        static constexpr lfsr_hash::STATE unit = {1, 0, 0, 0, 1, 0, 0, 0};
        const bool salt_ok = (g.g_251x4.get_state() != unit) &&
                             (g.g_241x4.get_state() != unit);

        report("v4: zero-salt does not reset state", power_ok && salt_ok,
               "state preserved", (power_ok && salt_ok) ? "state preserved" : "STATE WIPED");
    }

    /// Длина должна связываться полностью, а не только через чётность.
    void test_length_binding()
    {
        std::vector<uint8_t> a(64, 0x41), b(65, 0x41), zero(1, 0x00);
        const bool ok =
            lfsr_hash::hash128(std::as_bytes(std::span(a))) !=
                lfsr_hash::hash128(std::as_bytes(std::span(b))) &&
            lfsr_hash::hash128(std::span<const std::byte>{}) !=
                lfsr_hash::hash128(std::as_bytes(std::span(zero)));
        report("v4: full length is bound", ok, "different", ok ? "different" : "SAME");
    }
};
