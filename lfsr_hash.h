#pragma once

#include "lfsr.h"
#include <span>
#include <utility>

namespace lfsr_hash
{
    using LFSR251x4 = lfsr8::LFSR_paired_2x4<251>;
    using LFSR241x4 = lfsr8::LFSR_paired_2x4<241>;
    using STATE = lfsr8::u16x8;
    using u16 = lfsr8::u16;
    using u32 = lfsr8::u32;
    using u64 = lfsr8::u64;
    using u128 = std::pair<lfsr8::u64, lfsr8::u64>;

    static constexpr STATE K1 = {3, 2, 2, 0, 251-6, 5, 251-3, 5}; // p = 251.
    static constexpr STATE K2 = {7, 0, 4, 3, 241-7, 5, 1, 2};     // p = 241.

    struct salt
    {
        int q; // Количество тактов, применяемое к спаренному генератору.
        u16 s0; // Константный вход первого генератора.
        u16 s1; // -//- второго генератора.
    };

    static constexpr salt S0{7, 2, 3};
    static constexpr salt S1{6, 4, 7};
    static constexpr salt S2{31, 8, 11};
    static constexpr salt S3{29, 9, 5};
    static constexpr salt S4{37, 2, 13};

    static_assert(S0.q >= 4); // Минимально достаточное насыщение, m = 4
    static_assert(S1.q >= 4);

    static_assert(S2.q >= 6 * 4); // Насыщение с "запасом", m = 4
    static_assert(S3.q >= 6 * 4);
    static_assert(S4.q >= 6 * 4);

    struct gens
    {
        LFSR251x4 g_251x4; // Два спаренных генератора с периодами T1 = p1^4 - 1, T2 = p1^3 - 1: НОД(T1, T2) = p1 - 1.
        LFSR241x4 g_241x4; // Два спаренных генератора с периодами T1 = p2^4 - 1, T2 = p2^3 - 1: НОД(T1, T2) = p2 - 1.
        // Простые числа p1, p2 подобраны так, чтобы минимизировать НОД, но и стараться не выходить за пределы байта.

    public:
        /**
         * @brief Конструктор по умолчанию.
         */
        constexpr gens() : g_251x4(K1),
                           g_241x4(K2) {}
        
        /**
        * @brief Сбросить генераторы.
        * @details Устанавливается единичное состояние.
        */
        void reset()
        {
            g_251x4.set_unit_state();
            g_241x4.set_unit_state();
        }

        /**
         * @brief Добавить "соль". Проще говоря, по быстрому перевести генераторы в некоторое другое состояние,
         * не коррелированное с текущим.
         * @details Используется быстрое возведение в степень, т.к. кол-во тактов может быть большим.
         */
        inline void add_salt(const salt& S)
        {
            g_251x4.next_simd(S.s0);
            g_241x4.next_simd(S.s1);
            g_251x4.power_by(S.q);
            g_241x4.power_by(S.q);
        }

        /**
         * @brief Обработать входные байты с помощью генераторов.
         */
        void process_input(std::span<const std::byte> input);

        /**
         * @brief Сформировать 32-битный хеш исходя из текущих состояний генераторов.
         * @details Используется XOR состояний генераторов с некоторой (пока что фиксированной) маской.
         * Два спаренных генератора дают 4 состояния [x1, x2, x3, x4], размер которых 32 бита.
         */
        inline auto form_hash32()
        {
            auto st1 = g_251x4.get_state();
            auto st2 = g_241x4.get_state();
            lfsr8::u32 hash;
            hash = (st1[0] ^ 1 ^ st1[4]) ^ (st2[0] ^ 3 ^ st2[4]);
            hash <<= 8;
            hash |= (st1[1] ^ 3 ^ st1[5]) ^ (st2[1] ^ 1 ^ st2[5]);
            hash <<= 8;
            hash |= (st1[2] ^ 5 ^ st1[6]) ^ (st2[2] ^ 5 ^ st2[6]);
            hash <<= 8;
            hash |= (st1[3] ^ 3 ^ st1[7]) ^ (st2[3] ^ 3 ^ st2[7]);

            return hash;
        }
    };

    inline void lfsr_hash::gens::process_input(std::span<const std::byte> input)
    {
        const size_t n = input.size();
        if (n == 0) return;
        size_t i = 0;
        size_t j = n;

        for (; i + 16 <= j; i += 16)
        {
            // Читаем 16 байт целиком.
            __m128i data = _mm_loadu_si128((const __m128i *)(input.data() + i));

            g_251x4.next_simd_block(data);
            // Для второго генератора можно использовать инвертированный блок.
            g_241x4.next_simd_block(_mm_xor_si128(data, _mm_set1_epi8(0xFF)));
        }

        while (i < j)
        {
            uint16_t val_l = static_cast<uint8_t>(input[i]);
            uint16_t val_r = (i + 1 < j) ? static_cast<uint8_t>(input[j - 1]) : 0;

            g_251x4.next_simd(val_l, val_r);
            g_241x4.next_simd(~val_r, ~val_l);

            i++;
            if (j > i)
                j--;
        }

        uint8_t first = static_cast<uint8_t>(input[0]);
        uint16_t x = static_cast<uint16_t>(first | (first << 8));
        g_251x4.next_simd(x, x);
        g_241x4.next_simd(x, x);
    }

    /**
     * @brief Получить 32-битный хеш для входных байтов.
     */
    inline u32 hash32(gens &g, std::span<const std::byte> input)
    {
        const auto n = input.size();
        g.add_salt( n % 2 ? S1 : S0 );
        g.add_salt( n % 2 ? S0 : S1 );
        g.process_input(input);
        g.add_salt(S2);
        g.add_salt(S3);
        u32 h = g.form_hash32();
        return h;
    }

    /**
     * @brief Получить 64-битный хеш для входных байтов.
     */
    inline u64 hash64(gens &g, std::span<const std::byte> input)
    {
        const auto n = input.size();
        g.add_salt( n % 2 ? S1 : S0 );
        g.add_salt( n % 2 ? S0 : S1 );
        g.process_input(input);
        g.add_salt(S2);
        g.add_salt(S3);
        u64 h1 = g.form_hash32();
        g.add_salt(S2); // Нужна достаточная соль, чтобы удалиться от текущего состояния.
        g.add_salt(S4);
        u64 h2 = g.form_hash32();
        return (h1 << 32) | h2;
    }

    /**
     * @brief Получить 128-битный хеш для входных байтов.
     */
    inline u128 hash128(gens &g, std::span<const std::byte> input)
    {
        const auto n = input.size();
        g.add_salt( n % 2 ? S1 : S0 );
        g.add_salt( n % 2 ? S0 : S1 );
        g.process_input(input);
        g.add_salt(S1);
        g.add_salt(S3);

        u64 h1 = g.form_hash32();
        g.add_salt(S3); // Нужна достаточная соль, чтобы удалиться от текущего состояния.
        g.add_salt(S4);
        u64 h2 = g.form_hash32();
        g.add_salt(S2); // Нужна достаточная соль, чтобы удалиться от текущего состояния.
        g.add_salt(S4);
        u64 h3 = g.form_hash32();
        g.add_salt(S4); // Нужна достаточная соль, чтобы удалиться от текущего состояния.
        g.add_salt(S3);
        u64 h4 = g.form_hash32();
        return {
            h1 | (h2 << 32),
            h3 | (h4 << 32)};
    }
}