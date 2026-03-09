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
        const u32 diff = (state[0] >= input ? state[0] - input : state[0] - input + (SAMPLE)p) % (u16)p;
        const u32 m_v = ((u32)K_inv[0] * diff) % (u16)p;
        // 2. В цикле убираем лишние Modulo
        for (int i = 0; i < m - 1; i++)
        {
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
            if (i > 0 && baseT > ((std::numeric_limits<T>::max)() / result))
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
        void mult_by(const STATE &other)
        {
            STATE old_state = m_state;
            // Если не делать условную ссылку, то mult_by(m_state) даст неправильный результат, потому что
            // next(v) меняет m_state.
            const STATE &other_ref = other == m_state ? old_state : other;
            m_state.fill(0);
            for (int power = 2 * m - 2; power >= 0; --power)
            {
                u32 v = 0;
                for (int i = 0; i < power + 1; ++i)
                {
                    const int j = power - i;
                    if ((j >= m) || (j < 0))
                        continue;
                    if ((i >= m) || (i < 0))
                        continue;
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
        static STATE power_by(const STATE &state, const STATE &K, u64 q)
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
        static STATE inverse_of(const STATE &state, const STATE &K)
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
        bool is_state(const STATE &state) const
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
        constexpr LFSR_paired_2x4(u16x8 K) : m_K(K)
        {
            m_calculate_inverse_of_K();
            precompute_matrices();
        };

        void set_state(const u16x8 &state)
        {
            m_state = state;
        }

        void set_unit_state()
        {
            m_state = {1, 0, 0, 0, 1, 0, 0, 0};
        }

        void set_K(const u16x8 &K)
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

#if defined(_MSC_VER)
        __forceinline
#else
        [[gnu::always_inline]] inline
#endif
            void next_simd(u16 inp1, u16 inp2)
        {
            __m128i state = _mm_load_si128((const __m128i *)m_state.data());
            __m128i K_vec = _mm_load_si128((const __m128i *)m_K.data());

            // 1. Сдвиг влево (state[i] = state[i-1]):
            // В Little-endian памяти сдвиг к старшим индексам — это slli_si128
            __m128i shifted = _mm_slli_si128(state, 2);

            // 2. Изоляция: зануляем 0 и 4
            const __m128i iso_mask = _mm_setr_epi8(
                0, 0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                0, 0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff);
            shifted = _mm_and_si128(shifted, iso_mask);

            // 3. Feedback: v3 и v7. Используем _mm_set1_epi16 для правильного заполнения
            __m128i v_vec = _mm_set_epi16(m_state[7], m_state[7], m_state[7], m_state[7],
                                          m_state[3], m_state[3], m_state[3], m_state[3]);

            // 4. Входные данные: ставим inp1 в ячейку 0, inp2 в ячейку 4
            __m128i input_vec = _mm_setr_epi16(inp1, 0, 0, 0, inp2, 0, 0, 0);

            // 5. Расчет
            __m128i prod = _mm_mullo_epi16(v_vec, K_vec);
            __m128i res = _mm_add_epi16(shifted, _mm_add_epi16(prod, input_vec));

            _mm_store_si128((__m128i *)m_state.data(), simd_mod_p(res));
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

#if defined(_MSC_VER)
        __forceinline
#else
        [[gnu::always_inline]] inline
#endif
        void back_simd(u16 inp1, u16 inp2)
        {
            const u32 up = p;
            // 1. Восстанавливаем v1, v2 (обратная связь)
            // Важно: (state[0] - inp1 + p) % p — гарантируем отсутствие отрицательных чисел
            const u16 v1 = static_cast<u16>((u32)m_inv_K[0] * (m_state[0] + up - (inp1 % up)) % up);
            const u16 v2 = static_cast<u16>((u32)m_inv_K[4] * (m_state[4] + up - (inp2 % up)) % up);

            __m128i v_vec = _mm_setr_epi16(v1, v1, v1, v1, v2, v2, v2, v2);

            // 2. Загружаем K и сдвигаем его: нам нужны K[1], K[2], K[3] и K[5], K[6], K[7]
            // В цикле i=0..2: state[i] использует K[i+1]
            __m128i K_vec = _mm_load_si128((const __m128i *)m_K.data());
            __m128i K_plus_1 = _mm_srli_si128(K_vec, 2); // Сдвиг вправо на 1 элемент (2 байта)

            // 3. Загружаем state и сдвигаем: state[i+1]
            __m128i state = _mm_load_si128((const __m128i *)m_state.data());
            __m128i state_plus_1 = _mm_srli_si128(state, 2);

            // 4. Вычисление: (state[i+1] - v * K[i+1] + p*p) % p
            // Так как p < 256, v*K < 65536, mullo_epi16 здесь безопасен (без знакового переполнения)
            __m128i prod = _mm_mullo_epi16(v_vec, K_plus_1);
            __m128i prod_mod = simd_mod_p(prod);

            // (state[i+1] + p - prod_mod) % p
            __m128i p_vec = _mm_set1_epi16(static_cast<short>(p));
            __m128i diff = _mm_add_epi16(state_plus_1, _mm_sub_epi16(p_vec, prod_mod));
            __m128i res = simd_mod_p(diff);

            // 5. Маскируем и сохраняем
            // Нам нужно обновить только индексы 0,1,2 и 4,5,6.
            // Индексы 3 и 7 мы запишем отдельно.
            _mm_store_si128((__m128i *)m_state.data(), res);

            // 6. Финальная вставка v1 и v2 (те самые "новые" значения, пришедшие с конца)
            m_state[3] = v1;
            m_state[7] = v2;
        }

        // void next_simd_block(__m128i input_4_symbols_each) {
        //     // input_4_symbols_each должен содержать:
        //     // [inp1_t0, inp1_t1, inp1_t2, inp1_t3,  inp2_t0, inp2_t1, inp2_t2, inp2_t3]

        //     __m128i state = _mm_load_si128((const __m128i *)m_state.data());

        //     // 1. Предварительно вычисляем "вклады" каждой старой ячейки в новые.
        //     // За 4 шага s0, s1, s2, s3 полностью вытесняются и становятся "feedback"
        //     // для формирования новых s0...s3.

        //     // В векторизованном виде это эффективно делается через
        //     // умножение вектора состояния на заранее вычисленную матрицу переходов.

        //     // Однако, если p < 256, проще всего сделать 4 итерации внутри,
        //     // но используя SIMD для параллельной обработки ОБОИХ генераторов сразу.

        //     alignas(16) u16 inps[8];
        //     _mm_store_si128((__m128i*)inps, input_4_symbols_each);

        //     for(int t = 0; t < 4; ++t) {
        //         // Вызываем вашу исправленную next_simd,
        //         // передавая inps[t] (для 1-го) и inps[t+4] (для 2-го)
        //         next_simd(inps[t], inps[t+4]);
        //     }
        // }

        void next_simd_block(__m128i input_128)
        {
            __m128i s = _mm_load_si128((const __m128i *)m_state.data());

            // Бродкаст ячеек состояния (s0..s3 и s4..s7)
            __m128i bS[4];
            bS[0] = _mm_shuffle_epi8(s, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
            bS[1] = _mm_shuffle_epi8(s, _mm_setr_epi8(2, 3, 2, 3, 2, 3, 2, 3, 10, 11, 10, 11, 10, 11, 10, 11));
            bS[2] = _mm_shuffle_epi8(s, _mm_setr_epi8(4, 5, 4, 5, 4, 5, 4, 5, 12, 13, 12, 13, 12, 13, 12, 13));
            bS[3] = _mm_shuffle_epi8(s, _mm_setr_epi8(6, 7, 6, 7, 6, 7, 6, 7, 14, 15, 14, 15, 14, 15, 14, 15));

            // Бродкаст входов (inp_t0..t3)
            __m128i bI[4];
            bI[0] = _mm_shuffle_epi8(input_128, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 8, 9, 8, 9, 8, 9));
            bI[1] = _mm_shuffle_epi8(input_128, _mm_setr_epi8(2, 3, 2, 3, 2, 3, 2, 3, 10, 11, 10, 11, 10, 11, 10, 11));
            bI[2] = _mm_shuffle_epi8(input_128, _mm_setr_epi8(4, 5, 4, 5, 4, 5, 4, 5, 12, 13, 12, 13, 12, 13, 12, 13));
            bI[3] = _mm_shuffle_epi8(input_128, _mm_setr_epi8(6, 7, 6, 7, 6, 7, 6, 7, 14, 15, 14, 15, 14, 15, 14, 15));

            // Накопление результата: S' = sum(Si * Mi) + sum(Ii * Gi)
            // Т.к. p < 256, сумма 8-ми произведений (8 * 255 * 255 = 520,200) НЕ влезет в u16 (65535).
            // Поэтому делаем промежуточный модуль после 4-х сложений.

            __m128i acc = _mm_mullo_epi16(bS[0], m_matM[0]);
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bS[1], m_matM[1]));
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bS[2], m_matM[2]));
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bS[3], m_matM[3]));
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться

            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[0], m_matG[0]));
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[1], m_matG[1]));
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[2], m_matG[2]));
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[3], m_matG[3]));

            _mm_store_si128((__m128i *)m_state.data(), simd_mod_p(acc));
        }

        // /**
        //  * @brief Шаг назад на 4 такта сразу.
        //  * @param input_block_128 Вектор, содержащий входные символы для 4-х шагов.
        //  * Формат: [inp1_t0, inp1_t1, inp1_t2, inp1_t3,  inp2_t0, inp2_t1, inp2_t2, inp2_t3]
        //  */
        // void back_simd_block(__m128i input_block_128)
        // {
        //     // 1. Выгружаем входные данные во временный буфер для доступа в цикле
        //     alignas(16) u16 inps[8];
        //     _mm_store_si128(reinterpret_cast<__m128i*>(inps), input_block_128);

        //     // 2. Выполняем 4 шага назад.
        //     // Важно: идем от последнего входа к первому (t3 -> t2 -> t1 -> t0)
        //     for (int t = 3; t >= 0; --t)
        //     {
        //         // Используем ранее отлаженную back_simd (или её внутреннюю логику)
        //         // Передаем входные символы для текущего (обратного) шага
        //         back_simd(inps[t], inps[t + 4]);
        //     }
        // }

#if defined(_MSC_VER)
__forceinline 
#else
[[gnu::always_inline]] inline 
#endif
void back_simd_block(__m128i input_128) {
    __m128i state = _mm_load_si128((const __m128i*)m_state.data());
    const __m128i K_vec = _mm_load_si128((const __m128i*)m_K.data());
    const __m128i p_vec = _mm_set1_epi16(static_cast<short>(p));
    const __m128i K_shifted = _mm_srli_si128(K_vec, 2);

    // Используем макрос или лямбду с шаблонным параметром для индекса
    auto perform_back_step = [&](__m128i cur_s, const int t_idx, const int t_idx_p4) -> __m128i {
        // Чтобы обойти ограничение C2057, используем switch или ручной выбор
        // Но в развернутом коде проще всего вытащить значения заранее или через if constexpr
        u16 inp1, inp2;
        
        // Ручной выбор индекса (компилятор оптимизирует это в одну инструкцию pextrw)
        switch(t_idx) {
            case 0: inp1 = _mm_extract_epi16(input_128, 0); inp2 = _mm_extract_epi16(input_128, 4); break;
            case 1: inp1 = _mm_extract_epi16(input_128, 1); inp2 = _mm_extract_epi16(input_128, 5); break;
            case 2: inp1 = _mm_extract_epi16(input_128, 2); inp2 = _mm_extract_epi16(input_128, 6); break;
            case 3: default: inp1 = _mm_extract_epi16(input_128, 3); inp2 = _mm_extract_epi16(input_128, 7); break;
        }

        u16 s0 = static_cast<u16>(_mm_extract_epi16(cur_s, 0));
        u16 s4 = static_cast<u16>(_mm_extract_epi16(cur_s, 4));

        u16 v1 = static_cast<u16>((u32)m_inv_K[0] * (s0 + p - (inp1 % p)) % p);
        u16 v2 = static_cast<u16>((u32)m_inv_K[4] * (s4 + p - (inp2 % p)) % p);

        __m128i v_vec = _mm_setr_epi16(v1, v1, v1, v1, v2, v2, v2, v2);
        __m128i next_vals = _mm_srli_si128(cur_s, 2); 
        __m128i prod = simd_mod_p(_mm_mullo_epi16(v_vec, K_shifted));
        __m128i res = simd_mod_p(_mm_add_epi16(next_vals, _mm_sub_epi16(p_vec, prod)));

        // Вместо blend можно использовать _mm_insert_epi16 (тоже требует константу)
        // Но blend с константой 0x88 работает отлично
        __m128i v_insert = _mm_setr_epi16(0, 0, 0, v1, 0, 0, 0, v2);
        return _mm_blend_epi16(res, v_insert, 0x88);
    };

    // Разворачиваем шаги вручную, передавая константы
    state = perform_back_step(state, 3, 7);
    state = perform_back_step(state, 2, 6);
    state = perform_back_step(state, 1, 5);
    state = perform_back_step(state, 0, 4);

    _mm_store_si128((__m128i*)m_state.data(), state);
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
        void mult_by(const u16x8 &other)
        {
            u16x8 old_state = m_state;
            // Если не делать условную ссылку, то mult_by(m_state) даст неправильный результат, потому что
            // next(v) меняет m_state.
            const auto &other_ref = other == m_state ? old_state : other;
            m_state.fill(0);
            for (int power = 2 * 4 - 2; power >= 0; --power)
            {
                u32 v1 = 0;
                u32 v2 = 0;
                for (int i = 0; i < power + 1; ++i)
                {
                    const int j = power - i;
                    if ((j >= 4) || (j < 0))
                        continue;
                    if ((i >= 4) || (i < 0))
                        continue;
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
            for (; x != 0;)
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
        bool is_state_low(const u16x8 &state) const
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
        bool is_state_high(const u16x8 &state) const
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

    private:
        // Матрицы перехода для 4-х шагов (каждая строка — это влияние s_i на все s'_0..3)
        alignas(16) __m128i m_matM[4]; // Для состояния s0..s3 и s4..s7
        alignas(16) __m128i m_matG[4]; // Для входов inp_t0..t3

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

        void precompute_matrices()
        {
            // Для каждой ячейки i от 0 до 3 (считаем один генератор, второй симметричен)
            for (int i = 0; i < 4; ++i)
            {
                // 1. Влияние состояния: ставим 1 в ячейку i, остальные 0, входы 0
                m_state.fill(0);
                m_state[i] = 1;     // Для первого LFSR
                m_state[i + 4] = 1; // Для второго LFSR

                // Прогоняем 4 чистых шага (скалярных)
                for (int t = 0; t < 4; ++t)
                {
                    // Временная имитация next без внешних входов
                    // (используйте вашу эталонную логику next(0, 0) здесь)
                    next(0, 0);
                }
                // Результат после 4 шагов и есть i-я строка матрицы M
                m_matM[i] = _mm_loadu_si128((const __m128i *)&m_state[0]);
            }

            // 2. Влияние входов: подаем 1 на шаге t, в остальное время 0
            for (int t = 0; t < 4; ++t)
            {
                m_state.fill(0);
                for (int step = 0; step < 4; ++step)
                {
                    u16 inp = (step == t) ? 1 : 0;
                    next(inp, inp);
                }
                m_matG[t] = _mm_loadu_si128((const __m128i *)&m_state[0]);
            }
        }

        // Вспомогательная функция для векторизованного (x % p)
        // Работает сразу с 8 значениями u16
        inline __m128i simd_mod_p(__m128i x)
        {
            if constexpr (p == 256)
                return x; // Если p=256, ничего делать не надо

            // Константа Барретта: m = floor(2^16 / p)
            const __m128i m = _mm_set1_epi16(static_cast<short>((1u << 16) / p));
            const __m128i p_vec = _mm_set1_epi16(static_cast<short>(p));

            // q = (x * m) >> 16
            __m128i q = _mm_mulhi_epu16(x, m);
            // res = x - q * p
            __m128i res = _mm_sub_epi16(x, _mm_mullo_epi16(q, p_vec));

            // Коррекция: if (res >= p) res -= p
            __m128i mask = _mm_cmpgt_epi16(res, _mm_set1_epi16(static_cast<short>(p - 1)));
            return _mm_sub_epi16(res, _mm_and_si128(mask, p_vec));
        }
    };
}