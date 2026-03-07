#pragma once

/**
 * @author Новиков А.В., nawww83@gmail.com.
 *
 * Генератор LFSR в поле GF(p^m) с числом ячеек m = [1, 8] и простым модулем p = [2, 256*256).
 * Реализованы:
 *  - Генератор общего назначения:
 *     - p = [2, 256) для m = [5, 8],
 *     - p = [2, 256*256) для m = [1, 4].
 *  - Сдвоенный генератор с фиксированным m = 4, p = [2, 256). Фактически, это два независимых генератора.
 */

#include <cstdint>
#include <cassert>
#include <array>
#include <type_traits>
#include <cmath>
#include <concepts>
#include <limits>

#if defined(__x86_64__) || defined(_M_X64)
#define USE_SSE
#endif

#ifdef USE_SSE
#include <immintrin.h>
#include <smmintrin.h>
#endif

namespace lfsr8
{

    using u64 = uint64_t;
    using u32 = uint32_t;
    using u16 = uint16_t;
    using u16x8 = std::array<u16, 8>;
    using u32x4 = std::array<u32, 4>;

    template <int m>
    class MType
    {
    public:
        typedef typename std::conditional<(m <= 4), u32, u16>::type SAMPLE;
        typedef typename std::conditional<(m <= 4), u32x4, u16x8>::type STATE;
    };

    template <int p, int m>
    inline typename lfsr8::MType<m>::SAMPLE square_of_p()
    {
        using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
        return SAMPLE(p) * SAMPLE(p);
    }

    template <int p, int m>
    inline void modulo_p(typename lfsr8::MType<m>::STATE &state)
    {
        using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
        for (int i = 0; i < m; i++)
        {
            state[i] %= (SAMPLE)p;
        }
    }

    // Вспомогательная функция для скалярного Next (без SSE)
    template <int p, int m>
    inline void scalar_next(typename lfsr8::MType<m>::STATE &state, const typename lfsr8::MType<m>::STATE &K, typename lfsr8::MType<m>::SAMPLE input)
    {
        using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
        const SAMPLE m_v = state[m - 1];
        for (int i = m - 1; i > 0; i--)
        {
            state[i] = (state[i - 1] + m_v * K[i]) % (SAMPLE)p;
        }
        state[0] = (input + m_v * K[0]) % (SAMPLE)p;
    }

    // Вспомогательная функция для скалярного Back() (без SSE)
    template <int p, int m>
    inline void scalar_back(typename lfsr8::MType<m>::STATE &state, const typename lfsr8::MType<m>::STATE &K, const typename lfsr8::MType<m>::STATE &K_inv, typename lfsr8::MType<m>::SAMPLE input)
    {
        using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
        // 1. Считаем m_v один раз (здесь деление неизбежно)
        const u32 diff = ( state[0] >= input ? state[0] - input : state[0] - input + (SAMPLE)p ) % (u16)p;
        const u32 m_v = ( (u32)K_inv[0] * diff ) % (u16)p;
        // 2. В цикле убираем лишние Modulo
        for (int i = 0; i < m - 1; i++) {
            // Явно приводим к типу SAMPLE, который завязан на шаблонный p
            // Это заставит компилятор использовать оптимизированную последовательность для константы
            const u32 prod = m_v * (u32)K[i + 1];
            const u32 subtrahend = prod % (u16)p; // p — параметр шаблона
            // Финальное вычитание через условие (без деления!)
            int32_t res = (int32_t)state[i + 1] - (int32_t)subtrahend;
            res += res < 0 ? (int32_t)p : (int32_t)0;
            state[i] = (SAMPLE)res;
        }
        state[m - 1] = (SAMPLE)m_v;

        // Оригинальный вариант.
        // SAMPLE m_v = (state[0] - input + (SAMPLE)p) % (SAMPLE)p; // Берем по модулю перед последующим умножением,
        // т.к. может быть переполнение (т.е. неправильный модуль).
        // m_v *= K_inv[0];
        // m_v %= (SAMPLE)p;
        // for (int i = 0; i < m - 1; i++) {
            // state[i] = (state[i + 1] - m_v * K[i + 1] + (SAMPLE)(p) * (SAMPLE)(p)) % (SAMPLE)p;
        // }
        // state[m - 1] = m_v;
    }

    template <std::unsigned_integral T>
    constexpr T safe_ipow(int base, int exp)
    {
        T baseT = base;
        T result = 1;
        for (int i = 0; i < exp; ++i)
        {
            // Проверка: не превысит ли следующее умножение предел типа T
            if ( i > 0 && baseT > ((std::numeric_limits<T>::max)() / result) )
            {
                return 0; // Возвращаем 0 как маркер переполнения
            }
            result *= baseT;
        }
        return result;
    }

    /**
     * @brief Генератор LFSR общего назначения в поле GF(p^m).
     * @details
     * p должно быть простым числом в интервале:
     *   [2, 256) для длин регистра (4, 8].
     *   [2, 256*256) для длин регистра [1, 4].
     * m - длина регистра, [1, 8].
     */
    template <int p, int m>
    class LFSR
    {
        using STATE = typename MType<m>::STATE;
        using SAMPLE = typename MType<m>::SAMPLE;

    public:
        /**
         * @brief Конструктор с параметром.
         * @param K Коэффициенты (-g[0], ..., -g[m-1]) порождающего полинома g(x) степени m в поле GF(p).
         * @details Коэффициент при x^m равен 1 и не используется в явном виде.
         */
        constexpr LFSR(STATE K) : m_K(K)
        {
            static_assert(m <= 8);
            static_assert(m > 0);
            if constexpr (m > 4)
            {
                static_assert(p < 256);
            }
            else
            {
                static_assert(p < 256 * 256);
            }
            static_assert(p > 1);
            m_calculate_inverse_of_K();
        };

         /**
         * @brief Максимальный период генератора.
         */
        static constexpr auto T_MAX = safe_ipow<u64>(p, m) - 1;

        // Проверяем, что переполнения не случилось (P_POW_M не стал 0)
        static_assert(T_MAX != 0,
                      "Ошибка: Период превышает возможности 64-битного целого числа!");

        /**
         * @brief Установить состояние.
         * @param state Состояние.
         */
        void set_state(STATE state)
        {
            m_state = std::move(state);
        }

        /**
         * @brief Установить единичное состояние.
         */
        void set_unit_state()
        {
            m_state = {1};
        }

        /**
         * Установить коэффициенты порождающего полинома g(x).
         */
        void set_K(STATE K)
        {
            m_K = K;
            m_calculate_inverse_of_K();
        }

        /**
         * @brief Сделать шаг вперед (осуществить один такт генератора).
         * @param input Входной символ (должен быть приведен по модулю p!), который подается на вход генератора.
         */
        void next(SAMPLE input = 0)
        {
#ifdef USE_SSE
            if constexpr (m > 4)
            {
                __m128i a = _mm_set1_epi16(m_state[m - 1]);
                __m128i b = _mm_load_si128((const __m128i *)&m_K[0]);
                __m128i tmp = _mm_mullo_epi16(a, b);
                __m128i c = _mm_load_si128((const __m128i *)&m_state[0]);
                c = _mm_slli_si128(c, 2);
                __m128i mask = _mm_slli_si128(_mm_set1_epi16(-1), 2);
                __m128i inp = _mm_andnot_si128(mask, _mm_set1_epi16(input));
                c = _mm_add_epi16(c, tmp);
                c = _mm_add_epi16(c, inp);
                _mm_store_si128((__m128i *)&m_state[0], c);
                modulo_p<p, m>(m_state);
            }
            else
            {
                scalar_next<p, m>(m_state, m_K, input);
            }
#else // Для процессора общего назначения.
            scalar_next<p, m>(m_state, m_K, input);
#endif
        }

        /**
         * @brief Сделать шаг назад. Является операцией, обратной к next(input).
         * @param input Входной символ (должен быть приведен по модулю p!), который подается на вход генератора.
         */
        void back(SAMPLE input = 0)
        {
            scalar_back<p, m>(m_state, m_K, m_inv_K, input);
        }

        /**
         * @brief Возвести в квадрат, то есть вычислить состояние (x^s)^2, где
         * x^s - текущее состояние (некоторая степень s вспомогательной переменной x).
         * Соответствует s итерациям next() - прямое вычисление (долго, если s - большое число).
         */
        void square()
        {
           mult_by(m_state);
        }

        /**
         * @brief Умножить текущее состояние x^s на некоторое другое состояние x^t.
         * Итоговое состояние генератора становится равным x^(s+t).
         * @param other Другое состояние x^t.
         */
        void mult_by(const STATE& other)
        {
            STATE old_state = m_state;
            // Если не делать условную ссылку, то mult_by(m_state) даст неправильный результат, потому что
            // next(v) меняет m_state.
            const STATE& other_ref = other == m_state ? old_state : other;
            m_state.fill(0);
            for (int power = 2 * m - 2; power >= 0; --power)
            {
                u32 v = 0;
                for (int i = 0; i < power + 1; ++i)
                {
                    const int j = power - i;
                    if ((j >= m) || (j < 0)) continue;
                    if ((i >= m) || (i < 0)) continue;
                    v += ((u32)old_state[i] * (u32)other_ref[j]) % (u16)p;
                }
                next(static_cast<u16>(v));
            }
        }

        /**
         * @brief Возвести текущее состояние генератора в степень.
         * @param q Показатель степени.
         */
        void power_by(u64 q)
        {
            auto x = q;
            LFSR<p, m> lfsr{m_K};
            lfsr.set_unit_state();
            for (; x != 0;)
            {
                if ((x & 1) == 1)
                    lfsr.mult_by(get_state());
                square();
                x /= 2;
            }
            set_state(lfsr.get_state());
        }

        /**
         * @brief Возвести состояние в степень.
         * @param state Состояние x^t.
         * @param K Коэффициенты порождающего полинома g(x).
         * @param q Показатель степени.
         * @return Состояние x^(qt).
         */
        static STATE power_by(const STATE& state, const STATE& K, u64 q)
        {
            LFSR<p, m> lfsr{K};
            lfsr.set_state(state);
            lfsr.power_by(q);
            return lfsr.get_state();
        }

        /**
         * @brief Вычислить обратное состояние.
         * @param state Состояние x^t.
         * @param K Коэффициенты порождающего полинома g(x).
         * @return Состояние, 1/x^t.
         */
        static STATE inverse_of(const STATE& state, const STATE& K)
        {
            return power_by(state, K, T_MAX - 1);
        }

        /**
         * @brief Насытить генератор.
         * @param q Количество тактов.
         */
        void saturate(int q = m)
        {
            for (int i = 0; i < q; ++i)
            {
                next();
            }
        }

        /**
         * @brief Является ли заданное состояние текущим состоянием генератора.
         * @param state Заданное состояние.
         * @return Да/нет.
         */
        bool is_state(const STATE& state) const
        {
#ifdef USE_SSE
            bool result = true;
            for (int i = 0; i < m; ++i)
            {
                result &= (m_state[i] == state[i]);
            }
            return result;
#else
            return (state == m_state);
#endif
        }

        STATE get_generator_coeffs() const
        {
            return m_K;
        }

        STATE get_inverse_of_coeffs() const
        {
            return m_inv_K;
        }

        STATE get_state() const
        {
            return m_state;
        }

        SAMPLE get_cell(int idx) const
        {
            return m_state[idx];
        }

    private:
        /**
         * @brief Состояние генератора.
         */
        alignas(16) STATE m_state{};

        /**
         * @brief Коэффициенты порождающего полинома g(x).
         */
        alignas(16) STATE m_K{};

        /**
         * @brief Коэффициенты регистра для выполнения шага назад.
         */
        alignas(16) STATE m_inv_K{};

        /**
         * @brief Вычисляется обратный (по умножению) коэффициент K[0].
         */
        void m_calculate_inverse_of_K()
        {
            const auto g0 = m_K[0];
            assert(g0 != 0);
            SAMPLE inverse = 1;
            const SAMPLE modulo = p;
            for (;; inverse++)
            {
                const auto product = (g0 * inverse) % modulo;
                if (product == 1)
                    break;
            }
            m_inv_K[0] = inverse;
        }
    };

    /**
     * @brief Класс сдвоенного LFSR генератора общей длиной m = 4*2.
     * @details
     * Хранит числа в 16-битных ячейках.
     * Цель данного генератора: оптимизировать использование основного класса (см. LFSR), если
     * требуется парная работа генераторов. Генераторы работают независимо, но в
     * одном 128-битном регистре.
     */
    template <int p>
    class LFSR_paired_2x4
    {
        static_assert(p < 256);
        static_assert(p > 1);

    public:
        /**
         * @brief Конструктор с параметром.
         * @param K Коэффициенты (-g1[0], -g1[1], -g1[2], -g1[3], -g2[0], -g2[1], -g2[2], -g2[3]) двух
         * порождающих полиномов степени 4 в поле GF(p^4).
         * @details Коэффициент при x^m равен 1 и не используется в явном виде.
         */
        constexpr LFSR_paired_2x4(u16x8 K) : m_K(K) { m_calculate_inverse_of_K(); };

        void set_state(const u16x8& state)
        {
            m_state = state;
        }

        void set_unit_state()
        {
            m_state = {1, 0, 0, 0, 1, 0, 0, 0};
        }

        void set_K(const u16x8& K)
        {
            m_K = K;
            m_calculate_inverse_of_K();
        }

        /**
         * @brief Сделать шаг вперед (один такт генератора).
         * @param input Входной символ (должен быть приведен по модулю p!), который одинаково
         * подается на оба генератора.
         */
        void next(u16 input = 0)
        {
            next(input, input);
        }

        /**
         * @brief Сделать шаг вперед (один такт генератора).
         * @param inp1 Входной символ (должен быть приведен по модулю p!) первого генератора.
         * @param inp2 Входной символ (должен быть приведен по модулю p!) второго генератора.
         */
        void next(u16 inp1, u16 inp2)
        {
            u16 m_v3 = m_state[3];
            u16 m_v7 = m_state[7];
            for (int i = 7; i > 4; i--)
            {
                m_state[i] = ((u32)m_state[i - 1] + (u32)m_v7 * (u32)m_K[i]) % (u16)p;
                m_state[i - 4] = ((u32)m_state[i - 1 - 4] + (u32)m_v3 * (u32)m_K[i - 4]) % (u16)p;
            }
            m_state[0] = (inp1 + (u32)m_v3 * (u32)m_K[0]) % (u16)p;
            m_state[4] = (inp2 + (u32)m_v7 * (u32)m_K[4]) % (u16)p;
        }

        /**
         * @brief Сделать шаг назад (один такт генератора). Обратно к next(inp1, inp2).
         * @param inp1 Входной символ (должен быть приведен по модулю p!) первого генератора.
         * @param inp2 Входной символ (должен быть приведен по модулю p!) второго генератора.
         */
        void back(u16 inp1, u16 inp2)
        {
            const u16 m_v_1 = ((u32)m_inv_K[0] * ((u32)m_state[0] - (u32)inp1 + (u32)p)) % (u16)p;
            const u16 m_v_2 = ((u32)m_inv_K[4] * ((u32)m_state[4] - (u32)inp2 + (u32)p)) % (u16)p;
            for (int i = 0; i < 3; i++)
            {
                m_state[i] = ((u32)m_state[i + 1] - (u32)m_v_1 * (u32)m_K[i + 1] + (u32)p * (u32)p) % (u16)p;
                m_state[i + 4] = ((u32)m_state[i + 5] - (u32)m_v_2 * (u32)m_K[i + 5] + (u32)p * (u32)p) % (u16)p;
            }
            m_state[3] = m_v_1;
            m_state[7] = m_v_2;
        }

        /**
         * @brief Возвести в квадрат, то есть вычислить состояние (x^s)^2, где
         * x^s - текущее состояние (некоторая степень s вспомогательной переменной x).
         * Соответствует s итерациям next() - прямое вычисление (долго, если s - большое число).
         */
        void square()
        {
           mult_by(m_state);
        }

        /**
         * @brief Умножить текущее состояние x^s на некоторое другое состояние x^t.
         * Итоговое состояние генератора становится равным x^(s+t).
         * @param other Другое состояние x^t.
         */
        void mult_by(const u16x8& other)
        {
            u16x8 old_state = m_state;
            // Если не делать условную ссылку, то mult_by(m_state) даст неправильный результат, потому что
            // next(v) меняет m_state.
            const auto& other_ref = other == m_state ? old_state : other;
            m_state.fill(0);
            for (int power = 2 * 4 - 2; power >= 0; --power)
            {
                u32 v1 = 0;
                u32 v2 = 0;
                for (int i = 0; i < power + 1; ++i)
                {
                    const int j = power - i;
                    if ((j >= 4) || (j < 0)) continue;
                    if ((i >= 4) || (i < 0)) continue;
                    v1 += ((u32)old_state[i] * (u32)other_ref[j]) % (u16)p;
                    v2 += ((u32)old_state[i + 4] * (u32)other_ref[j + 4]) % (u16)p;
                }
                next(static_cast<u16>(v1), static_cast<u16>(v2));
            }
        }

        void power_by(u64 q)
        {
            auto x = q;
            LFSR_paired_2x4<p> lfsr{m_K};
            lfsr.set_unit_state();
            for (; x != 0; )
            {
                if ((x & 1) == 1)
                    lfsr.mult_by(get_state());
                square();
                x /= 2;
            }
            set_state(lfsr.get_state());
        }

        auto get_state() const
        {
            return m_state;
        }

        /**
         * @brief Совпадает ли заданное состояние с текущей нижней частью состояния.
         * @param state Заданное состояние.
         * @return Да/нет.
         */
        bool is_state_low(const u16x8& state) const
        {
            bool result = true;
            for (int i = 0; i < 4; ++i)
            {
                result &= (m_state[i] == state[i]);
            }
            return result;
        }

        /**
         * @brief Совпадает ли заданное состояние с текущей верхней частью состояния.
         * @param state Заданное состояние.
         * @return Да/нет.
         */
        bool is_state_high(const u16x8& state) const
        {
            bool result = true;
            for (int i = 0; i < 4; ++i)
            {
                result &= (m_state[i + 4] == state[i + 4]);
            }
            return result;
        }

    private:
        alignas(16) u16x8 m_state{};

        alignas(16) u16x8 m_K{};

        alignas(16) u16x8 m_inv_K{};

        /**
         * @brief Вычисляются обратные (по умножению) коэффициенты.
         */
        void m_calculate_inverse_of_K()
        {
            const u32 g0 = m_K[0]; // Берем u32 для защиты от переполнения при умножении.
            const u32 g4 = m_K[4];
            assert(g0 != 0);
            assert(g4 != 0);
            u32 inverse0 = 1;
            u32 inverse4 = 1;
            const u32 modulo = p;
            bool achieved0 = false;
            bool achieved4 = false;
            for (;;)
            {
                const auto product0 = (g0 * inverse0) % modulo;
                const auto product4 = (g4 * inverse4) % modulo;
                achieved0 = achieved0 ? achieved0 : product0 == 1;
                achieved4 = achieved4 ? achieved4 : product4 == 1;
                if (achieved0 && achieved4)
                    break;
                inverse0 += !achieved0;
                inverse4 += !achieved4;
            }
            m_inv_K[0] = static_cast<u16>(inverse0);
            m_inv_K[4] = static_cast<u16>(inverse4);
        }
    };
}