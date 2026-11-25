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
#include <cstddef>
#include <cassert>
#include <array>
#include <type_traits>
#include <concepts>
#include <limits>

#if defined(__x86_64__) || defined(_M_X64)
#define USE_SSE
#endif

#ifdef USE_SSE
#include <immintrin.h>
#include <smmintrin.h>
#else
#error "lfsr.h: SIMD-ядро реализовано только для x86_64 (SSE4.2/AVX2); скалярный путь для других архитектур отсутствует."
#endif

// Слитный блоковый процессор (LFSR_fused_2x4) требует AVX2.
#if defined(USE_SSE) && defined(__AVX2__)
#define USE_AVX2
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
    inline void modulo_p(typename lfsr8::MType<m>::STATE &state)
    {
        using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
        for (int i = 0; i < m; i++)
        {
            state[i] %= (u32)p;
        }
    }

    // Вспомогательная функция для скалярного next()
    template <int p, int m>
    inline void scalar_next(typename lfsr8::MType<m>::STATE &state, const typename lfsr8::MType<m>::STATE &K, typename lfsr8::MType<m>::SAMPLE input)
    {
        using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
        const SAMPLE m_v = state[m - 1];
        for (int i = m - 1; i > 0; i--)
        {
            state[i] = (state[i - 1] + (u32)m_v * K[i]) % (u32)p;
        }
        state[0] = (input + (u32)m_v * K[0]) % (u32)p;
    }

    // Вспомогательная функция для скалярного back()
    template <int p, int m>
    inline void scalar_back(typename lfsr8::MType<m>::STATE &state, const typename lfsr8::MType<m>::STATE &K, const typename lfsr8::MType<m>::STATE &K_inv, typename lfsr8::MType<m>::SAMPLE input)
    {
        using SAMPLE = typename lfsr8::MType<m>::SAMPLE;
        // Приводим вход и текущую ячейку по модулю p, чтобы не зависеть от
        // переполнения беззнакового типа при input >= p или невыровненном состоянии.
        const u32 up = (u32)p;
        const u32 in = (u32)input % up;
        const u32 cur = (u32)state[0] % up;
        const u32 diff = cur >= in ? cur - in : cur + up - in;
        const u32 m_v = ((u32)K_inv[0] * diff) % up;
        for (int i = 0; i < m - 1; i++)
        {
            const u32 prod = m_v * (u32)K[i + 1];
            const u32 subtrahend = prod % (u32)p;
            int32_t res = (int32_t)state[i + 1] - (int32_t)subtrahend;
            res += res < 0 ? (int32_t)p : (int32_t)0;
            state[i] = (SAMPLE)res;
        }
        state[m - 1] = (SAMPLE)m_v;
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
         * @param input Входной символ, который подается на вход генератора (приводится по модулю p внутри).
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
         * @param input Входной символ, который подавался на вход генератора (приводится по модулю p внутри).
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
                    v += ((u32)old_state[i] * other_ref[j]) % (u32)p;
                }
                next(v % (u32)p);
            }
        }

        /**
         * @brief Возвести текущее состояние генератора в степень.
         * @param q Показатель степени.
         * @details q = 0 — тождество (состояние не меняется). Формально x^0 = 1,
         * но для генератора это означало бы схлопывание в единичное состояние,
         * одинаковое для всех сообщений, поэтому нулевая степень трактуется как
         * «сдвиг не применять».
         */
        void power_by(u64 q)
        {
            if (q == 0)
                return;
            auto x = q;
            LFSR<p, m> lfsr{m_K};
            lfsr.set_unit_state();
            for (;;)
            {
                if ((x & 1) == 1)
                    lfsr.mult_by(get_state());
                x >>= 1;
                if (x == 0)
                    break; // последний square() уже не нужен
                square();
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
         * @param input Входной символ, который одинаково подается на оба генератора.
         */
        void next(u16 input = 0)
        {
            next(input, input);
        }

        /**
         * @brief Сделать шаг вперед (один такт генератора).
         * @param inp1 Входной символ первого генератора.
         * @param inp2 Входной символ второго генератора.
         */
        void next(u16 inp1, u16 inp2)
        {
            const u32 up = p;
            u16 m_v3 = m_state[3];
            u16 m_v7 = m_state[7];
            for (int i = 7; i > 4; i--)
            {
                m_state[i - 0] = ((u32)m_state[i - 1 - 0] + (u32)m_v7 * m_K[i - 0]) % up;
                m_state[i - 4] = ((u32)m_state[i - 1 - 4] + (u32)m_v3 * m_K[i - 4]) % up;
            }
            m_state[0] = (inp1 + (u32)m_v3 * m_K[0]) % up;
            m_state[4] = (inp2 + (u32)m_v7 * m_K[4]) % up;
        }

        void next_simd(u16 inp)
        {
            next_simd(inp, inp);
        }

#if defined(_MSC_VER)
        __forceinline
#else
        [[gnu::always_inline]] inline
#endif
            void next_simd(u16 inp1, u16 inp2)
        {
            // Вход приводится по модулю p: иначе сумма state + feedback + input
            // может переполнить 16 бит (так происходит, например, для
            // дополненных байтов в хвосте) и результат разойдётся с next().
            inp1 = static_cast<u16>(inp1 % p);
            inp2 = static_cast<u16>(inp2 % p);

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

        void back(u16 inp)
        {
            back(inp, inp);
        }

        /**
         * @brief Сделать шаг назад (один такт генератора). Обратно к next(inp1, inp2).
         * @param inp1 Входной символ первого генератора.
         * @param inp2 Входной символ второго генератора.
         */
        void back(u16 inp1, u16 inp2)
        {
            const u32 up = p;
            const u16 m_v_1 = ((u32)m_inv_K[0] * ((u32)m_state[0] - (inp1 % up) + up)) % up;
            const u16 m_v_2 = ((u32)m_inv_K[4] * ((u32)m_state[4] - (inp2 % up) + up)) % up;
            for (int i = 0; i < 3; i++)
            {
                u32 a = ((u32)m_v_1 * m_K[i + 1]) % up;
                u32 b = ((u32)m_v_2 * m_K[i + 5]) % up;
                m_state[i + 0] = ((u32)m_state[i + 1] - a + up) % up;
                m_state[i + 4] = ((u32)m_state[i + 5] - b + up) % up;
            }
            m_state[3] = m_v_1;
            m_state[7] = m_v_2;
        }

        void back_simd(u16 inp)
        {
            back_simd(inp, inp);
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

        void next_simd_block(__m128i input_128)
        {
            __m128i s = _mm_load_si128((const __m128i *)m_state.data());

            // Входные слова приводятся по модулю p ДО умножения на матрицу G.
            // Иначе произведение input * G (до 65535 * 250) переполняет 16 бит
            // и результат расходится со скалярной моделью next().
            input_128 = simd_mod_p(input_128);

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
            __m128i acc = _mm_mullo_epi16(bS[0], m_matM[0]);
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bS[1], m_matM[1]));
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bS[2], m_matM[2]));
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bS[3], m_matM[3]));
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[0], m_matG[0]));
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[1], m_matG[1]));
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[2], m_matG[2]));
            acc = simd_mod_p(acc); // Промежуточный модуль, чтобы не переполниться
            acc = _mm_add_epi16(acc, _mm_mullo_epi16(bI[3], m_matG[3]));
            _mm_store_si128((__m128i *)m_state.data(), simd_mod_p(acc));
        }

#if defined(_MSC_VER)
        __forceinline
#else
        [[gnu::always_inline]] inline
#endif
            void back_simd_block(__m128i input_128)
        {
            __m128i state = _mm_load_si128((const __m128i *)m_state.data());
            const __m128i K_vec = _mm_load_si128((const __m128i *)m_K.data());
            const __m128i p_vec = _mm_set1_epi16(static_cast<short>(p));
            const __m128i K_shifted = _mm_srli_si128(K_vec, 2);

            // Используем макрос или лямбду с шаблонным параметром для индекса
            auto perform_back_step = [&](__m128i cur_s, const int t_idx, const int t_idx_p4) -> __m128i
            {
                // Чтобы обойти ограничение C2057, используем switch или ручной выбор
                // Но в развернутом коде проще всего вытащить значения заранее или через if constexpr
                u16 inp1, inp2;

                // Ручной выбор индекса (компилятор оптимизирует это в одну инструкцию pextrw)
                switch (t_idx)
                {
                case 0:
                    inp1 = _mm_extract_epi16(input_128, 0);
                    inp2 = _mm_extract_epi16(input_128, 4);
                    break;
                case 1:
                    inp1 = _mm_extract_epi16(input_128, 1);
                    inp2 = _mm_extract_epi16(input_128, 5);
                    break;
                case 2:
                    inp1 = _mm_extract_epi16(input_128, 2);
                    inp2 = _mm_extract_epi16(input_128, 6);
                    break;
                case 3:
                default:
                    inp1 = _mm_extract_epi16(input_128, 3);
                    inp2 = _mm_extract_epi16(input_128, 7);
                    break;
                }

                u16 s0 = static_cast<u16>(_mm_extract_epi16(cur_s, 0));
                u16 s4 = static_cast<u16>(_mm_extract_epi16(cur_s, 4));

                const u32 up = p;
                u16 v1 = static_cast<u16>((u32)m_inv_K[0] * (s0 + up - (inp1 % up)) % up);
                u16 v2 = static_cast<u16>((u32)m_inv_K[4] * (s4 + up - (inp2 % up)) % up);

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

            _mm_store_si128((__m128i *)m_state.data(), state);
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

        void mult_by(const u16x8 &other)
        {
            // 1. Делаем локальную копию ТЕКУЩЕГО состояния
            u16x8 a_copy = m_state;

            // 2. Решаем проблему square(): если нам передали ссылку на самого себя,
            // используем нашу копию a_copy. Если передали другой объект — используем его напрямую.
            const u16x8 &b_ref = (&other == &m_state) ? a_copy : other;

            // 3. Теперь можно безопасно обнулять m_state для накопления результата
            m_state.fill(0);

            for (int power = 6; power >= 0; --power)
            {
                uint32_t v1 = 0;
                uint32_t v2 = 0;

                int i_start = (std::max)(0, power - 3);
                int i_end = (std::min)(power, 3);

                for (int i = i_start; i <= i_end; ++i)
                {
                    // Используем a_copy (множимое) и b_ref (множитель)
                    v1 += (uint32_t)a_copy[i + 0] * b_ref[power - i + 0];
                    v2 += (uint32_t)a_copy[i + 4] * b_ref[power - i + 4];
                }

                // Выполняем шаг LFSR (редукция полинома)
                next_simd(static_cast<u16>(v1 % p),
                          static_cast<u16>(v2 % p));
            }
        }

        /**
         * @brief Возвести текущее состояние в степень q.
         * @details q = 0 — тождество (состояние не меняется). Формально x^0 = 1,
         * но для генератора это означало бы схлопывание в единичное состояние,
         * одинаковое для всех сообщений: соль с q = 0 схлопнула бы состояние
         * в единичное.
         */
        void power_by(u64 q)
        {
            if (q == 0)
                return;
            auto x = q;
            // Копия текущего объекта, чтобы не пересчитывать матрицы и
            // обратные коэффициенты в конструкторе на каждый вызов.
            LFSR_paired_2x4<p> lfsr{*this};
            lfsr.set_unit_state();
            for (;;)
            {
                if ((x & 1) == 1)
                    lfsr.mult_by(get_state());
                x >>= 1;
                if (x == 0)
                    break; // последний square() уже не нужен
                square();
            }
            set_state(lfsr.get_state());
        }

        auto get_state() const
        {
            return m_state;
        }

        /**
         * @brief Коэффициенты порождающего полинома g(x).
         */
        u16x8 get_generator_coeffs() const
        {
            return m_K;
        }

        /**
         * @brief Прямой доступ к ячейкам состояния.
         * @details Нужен слитному блоковому процессору, чтобы не копировать
         * состояние туда-обратно на каждый блок.
         */
        u16x8 &state_ref()
        {
            return m_state;
        }

        const u16x8 &state_ref() const
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
                m_state[i + 0] = 1; // Для первого LFSR
                m_state[i + 4] = 1; // Для второго LFSR

                // Прогоняем 4 чистых шага (скалярных)
                for (int t = 0; t < 4; ++t)
                    next(0, 0);
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

    /**
     * @brief Слитный блоковый процессор для двух спаренных генераторов с разными модулями.
     * @details
     * За один блок (BLOCK_BYTES = 512 байт) каждый генератор делает STEPS = 128
     * тактов. Вклады состояния и входов считаются через предвычисленные
     * матрицы перехода (A^STEPS и A^k), поэтому редукция по модулю выполняется
     * один раз на блок, а не на каждом такте.
     *
     * Матрицы строятся по скалярной модели next(), а входные слова приводятся
     * по модулю p, поэтому результат побитно совпадает с последовательным
     * скалярным вычислением (без 16-битного переполнения).
     *
     * Раскладка входа соответствует LFSR_paired_2x4::next_simd_block():
     * за 16 байт данных выполняется 4 такта, младшая половина слова идёт в
     * нижний 4-ячеечный регистр, старшая — в верхний; для второго генератора
     * байты предварительно инвертируются.
     *
     * AVX2: состояния обоих генераторов лежат в одном 256-битном регистре
     * (младшие 128 бит — p1, старшие — p2).
     * SSE4.2: те же самые матрицы, но каждый генератор считается в своём
     * отдельном 128-битном регистре. Математика идентична, поэтому обе сборки
     * дают побитно одинаковый хэш.
     */
    template <int p1, int p2>
    class LFSR_fused_2x4
    {
        static_assert(p1 > 1 && p1 < 256);
        static_assert(p2 > 1 && p2 < 256);

    public:
        /// Количество тактов генератора в одном блоке.
        static constexpr int STEPS = 128;

        /// Количество байт данных, обрабатываемых за один блок (4 байта на такт).
        static constexpr int BLOCK_BYTES = 4 * STEPS;

        constexpr LFSR_fused_2x4(const u16x8 &K1, const u16x8 &K2)
        {
            build(K1, K2);
        }

        /**
         * @brief Обработать nblocks блоков по BLOCK_BYTES байт.
         * @param data Входные данные.
         * @param nblocks Количество полных блоков.
         * @param state1 Состояние первого генератора (p1), обновляется на месте.
         * @param state2 Состояние второго генератора (p2), обновляется на месте.
         */
        void process_blocks(const std::byte *data, size_t nblocks,
                            u16x8 &state1, u16x8 &state2) const
        {
#ifdef USE_AVX2
            __m256i st = _mm256_set_m128i(
                _mm_load_si128(reinterpret_cast<const __m128i *>(state2.data())),
                _mm_load_si128(reinterpret_cast<const __m128i *>(state1.data())));

            for (size_t b = 0; b < nblocks; ++b)
            {
                const std::byte *src = data + b * BLOCK_BYTES;
                __m256i acc_lo = _mm256_setzero_si256();
                __m256i acc_hi = _mm256_setzero_si256();

                for (int k = 0; k < STEPS / 4; ++k)
                {
                    __m128i blk = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + 16 * k));
                    // Старшая половина инвертируется: второй генератор
                    // обрабатывает дополнение входных байтов.
                    __m256i x = _mm256_xor_si256(_mm256_broadcastsi128_si256(blk), m_mats.hi_inv);
                    x = reduce_p16(x);

                    const __m256i i01 = _mm256_shuffle_epi8(x, m_mats.m01);
                    const __m256i i23 = _mm256_shuffle_epi8(x, m_mats.m23);
                    const __m256i i45 = _mm256_shuffle_epi8(x, m_mats.m45);
                    const __m256i i67 = _mm256_shuffle_epi8(x, m_mats.m67);

                    acc_lo = _mm256_add_epi32(acc_lo, _mm256_madd_epi16(i01, m_mats.Glo[2 * k]));
                    acc_lo = _mm256_add_epi32(acc_lo, _mm256_madd_epi16(i23, m_mats.Glo[2 * k + 1]));
                    acc_hi = _mm256_add_epi32(acc_hi, _mm256_madd_epi16(i45, m_mats.Ghi[2 * k]));
                    acc_hi = _mm256_add_epi32(acc_hi, _mm256_madd_epi16(i67, m_mats.Ghi[2 * k + 1]));
                }

                const __m256i s01 = _mm256_shuffle_epi8(st, m_mats.m01);
                const __m256i s23 = _mm256_shuffle_epi8(st, m_mats.m23);
                const __m256i s45 = _mm256_shuffle_epi8(st, m_mats.m45);
                const __m256i s67 = _mm256_shuffle_epi8(st, m_mats.m67);

                const __m256i s_lo = _mm256_add_epi32(acc_lo,
                    _mm256_add_epi32(_mm256_madd_epi16(s01, m_mats.A01lo),
                                     _mm256_madd_epi16(s23, m_mats.A23lo)));
                const __m256i s_hi = _mm256_add_epi32(acc_hi,
                    _mm256_add_epi32(_mm256_madd_epi16(s45, m_mats.A01hi),
                                     _mm256_madd_epi16(s67, m_mats.A23hi)));

                st = _mm256_packus_epi32(mod_p32(s_lo), mod_p32(s_hi));
            }

            _mm_store_si128(reinterpret_cast<__m128i *>(state1.data()), _mm256_castsi256_si128(st));
            _mm_store_si128(reinterpret_cast<__m128i *>(state2.data()), _mm256_extracti128_si256(st, 1));
#else
            // SSE4.2: тот же алгоритм, каждый генератор — в своём регистре.
            __m128i st1 = _mm_load_si128(reinterpret_cast<const __m128i *>(state1.data()));
            __m128i st2 = _mm_load_si128(reinterpret_cast<const __m128i *>(state2.data()));
            const __m128i zero = _mm_setzero_si128();

            for (size_t b = 0; b < nblocks; ++b)
            {
                const std::byte *src = data + b * BLOCK_BYTES;
                __m128i a1lo = zero, a1hi = zero, a2lo = zero, a2hi = zero;

                for (int k = 0; k < STEPS / 4; ++k)
                {
                    const __m128i blk = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + 16 * k));
                    const __m128i x1 = reduce_p16(blk, m_p1);
                    const __m128i x2 = reduce_p16(_mm_xor_si128(blk, m_p2.ff), m_p2);

                    a1lo = _mm_add_epi32(a1lo, _mm_madd_epi16(_mm_shuffle_epi8(x1, m_p1.m01), m_p1.Glo[2 * k]));
                    a1lo = _mm_add_epi32(a1lo, _mm_madd_epi16(_mm_shuffle_epi8(x1, m_p1.m23), m_p1.Glo[2 * k + 1]));
                    a1hi = _mm_add_epi32(a1hi, _mm_madd_epi16(_mm_shuffle_epi8(x1, m_p1.m45), m_p1.Ghi[2 * k]));
                    a1hi = _mm_add_epi32(a1hi, _mm_madd_epi16(_mm_shuffle_epi8(x1, m_p1.m67), m_p1.Ghi[2 * k + 1]));

                    a2lo = _mm_add_epi32(a2lo, _mm_madd_epi16(_mm_shuffle_epi8(x2, m_p2.m01), m_p2.Glo[2 * k]));
                    a2lo = _mm_add_epi32(a2lo, _mm_madd_epi16(_mm_shuffle_epi8(x2, m_p2.m23), m_p2.Glo[2 * k + 1]));
                    a2hi = _mm_add_epi32(a2hi, _mm_madd_epi16(_mm_shuffle_epi8(x2, m_p2.m45), m_p2.Ghi[2 * k]));
                    a2hi = _mm_add_epi32(a2hi, _mm_madd_epi16(_mm_shuffle_epi8(x2, m_p2.m67), m_p2.Ghi[2 * k + 1]));
                }

                a1lo = _mm_add_epi32(a1lo, _mm_madd_epi16(_mm_shuffle_epi8(st1, m_p1.m01), m_p1.A01lo));
                a1lo = _mm_add_epi32(a1lo, _mm_madd_epi16(_mm_shuffle_epi8(st1, m_p1.m23), m_p1.A23lo));
                a1hi = _mm_add_epi32(a1hi, _mm_madd_epi16(_mm_shuffle_epi8(st1, m_p1.m45), m_p1.A01hi));
                a1hi = _mm_add_epi32(a1hi, _mm_madd_epi16(_mm_shuffle_epi8(st1, m_p1.m67), m_p1.A23hi));

                a2lo = _mm_add_epi32(a2lo, _mm_madd_epi16(_mm_shuffle_epi8(st2, m_p2.m01), m_p2.A01lo));
                a2lo = _mm_add_epi32(a2lo, _mm_madd_epi16(_mm_shuffle_epi8(st2, m_p2.m23), m_p2.A23lo));
                a2hi = _mm_add_epi32(a2hi, _mm_madd_epi16(_mm_shuffle_epi8(st2, m_p2.m45), m_p2.A01hi));
                a2hi = _mm_add_epi32(a2hi, _mm_madd_epi16(_mm_shuffle_epi8(st2, m_p2.m67), m_p2.A23hi));

                st1 = _mm_packus_epi32(mod_p32(a1lo, m_p1), mod_p32(a1hi, m_p1));
                st2 = _mm_packus_epi32(mod_p32(a2lo, m_p2), mod_p32(a2hi, m_p2));
            }

            _mm_store_si128(reinterpret_cast<__m128i *>(state1.data()), st1);
            _mm_store_si128(reinterpret_cast<__m128i *>(state2.data()), st2);
#endif
        }

    private:
        /// Вклады ячеек и входов за STEPS тактов (по одному генератору).
        struct trans
        {
            u16 Alo[4][4]{};        // Влияние ячеек нижнего регистра.
            u16 Ahi[4][4]{};        // Влияние ячеек верхнего регистра.
            u16 Glo[STEPS][4]{};    // Влияние входов нижнего регистра.
            u16 Ghi[STEPS][4]{};    // Влияние входов верхнего регистра.
        };

#ifdef USE_AVX2
        /// Предвычисленные векторы для слитного прохода.
        struct mats
        {
            __m256i A01lo, A23lo, A01hi, A23hi;
            __m256i Glo[STEPS / 2], Ghi[STEPS / 2];
            __m256i P, Pm1, R32, R32odd; // 32-битная редукция по модулю.
            __m256i P16, M16;            // 16-битная редукция входа.
            __m256i m01, m23, m45, m67;  // pshufb: пара соседних ячеек.
            __m256i hi_inv;              // Инверсия старшей половины.
        };

        alignas(32) mats m_mats{};
#endif

        /// Матрицы перехода одного генератора с модулем p.
        /// @details Строятся как степени матрицы перехода A (4x4), а не
        /// прогоном STEPS шагов для каждой ячейки и каждого входа: это O(STEPS)
        /// вместо O(STEPS^2) и заметно дешевле при больших блоках.
        ///   s'[0] = input + K[0]*s[3]
        ///   s'[i] = s[i-1] + K[i]*s[3], i = 1..3
        template <int p>
        static trans build_trans(const u16x8 &K)
        {
            trans r{};
            const u32 up = static_cast<u32>(p);

            for (int half = 0; half < 2; ++half)
            {
                const u16 *k = &K[half * 4];

                // B = A^step, v = A^step * e0.
                u32 B[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
                u32 v[4] = {1, 0, 0, 0};

                for (int step = 0; step < STEPS; ++step)
                {
                    // Вход на шаге t даёт вклад A^(STEPS-1-t) * e0.
                    const int t = STEPS - 1 - step;
                    for (int j = 0; j < 4; ++j)
                    {
                        if (half == 0)
                            r.Glo[t][j] = static_cast<u16>(v[j]);
                        else
                            r.Ghi[t][j] = static_cast<u16>(v[j]);
                    }

                    u32 nv[4] = {};
                    nv[0] = (static_cast<u32>(k[0]) * v[3]) % up;
                    for (int i = 1; i < 4; ++i)
                        nv[i] = (v[i - 1] + static_cast<u32>(k[i]) * v[3]) % up;
                    for (int i = 0; i < 4; ++i)
                        v[i] = nv[i];

                    u32 nb[4][4] = {};
                    for (int j = 0; j < 4; ++j)
                    {
                        nb[0][j] = (static_cast<u32>(k[0]) * B[3][j]) % up;
                        for (int i = 1; i < 4; ++i)
                            nb[i][j] = (B[i - 1][j] + static_cast<u32>(k[i]) * B[3][j]) % up;
                    }
                    for (int i = 0; i < 4; ++i)
                        for (int j = 0; j < 4; ++j)
                            B[i][j] = nb[i][j];
                }

                // Alo[i][j] — вклад исходной ячейки i в итоговую ячейку j,
                // то есть B[j][i] (B = A^STEPS).
                for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                    {
                        if (half == 0)
                            r.Alo[i][j] = static_cast<u16>(B[j][i]);
                        else
                            r.Ahi[i][j] = static_cast<u16>(B[j][i]);
                    }
            }
            return r;
        }

#ifdef USE_AVX2
        void build(const u16x8 &K1, const u16x8 &K2)
        {
            const trans t1 = build_trans<p1>(K1);
            const trans t2 = build_trans<p2>(K2);

            alignas(32) u16 buf[16];

            // Вектор из пар (A[base][j], A[base + 1][j]): дорожки 0..7 — p1,
            // дорожки 8..15 — p2.
            auto mk_pair = [&](auto getA, int base)
            {
                for (int j = 0; j < 4; ++j)
                {
                    buf[j * 2] = getA(t1, base, j);
                    buf[j * 2 + 1] = getA(t1, base + 1, j);
                    buf[8 + j * 2] = getA(t2, base, j);
                    buf[8 + j * 2 + 1] = getA(t2, base + 1, j);
                }
                return _mm256_load_si256(reinterpret_cast<const __m256i *>(buf));
            };

            m_mats.A01lo = mk_pair([](const trans &t, int i, int j) { return t.Alo[i][j]; }, 0);
            m_mats.A23lo = mk_pair([](const trans &t, int i, int j) { return t.Alo[i][j]; }, 2);
            m_mats.A01hi = mk_pair([](const trans &t, int i, int j) { return t.Ahi[i][j]; }, 0);
            m_mats.A23hi = mk_pair([](const trans &t, int i, int j) { return t.Ahi[i][j]; }, 2);

            for (int k = 0; k < STEPS / 2; ++k)
            {
                const int t = 2 * k;
                for (int j = 0; j < 4; ++j)
                {
                    buf[j * 2] = t1.Glo[t][j];
                    buf[j * 2 + 1] = t1.Glo[t + 1][j];
                    buf[8 + j * 2] = t2.Glo[t][j];
                    buf[8 + j * 2 + 1] = t2.Glo[t + 1][j];
                }
                m_mats.Glo[k] = _mm256_load_si256(reinterpret_cast<const __m256i *>(buf));

                for (int j = 0; j < 4; ++j)
                {
                    buf[j * 2] = t1.Ghi[t][j];
                    buf[j * 2 + 1] = t1.Ghi[t + 1][j];
                    buf[8 + j * 2] = t2.Ghi[t][j];
                    buf[8 + j * 2 + 1] = t2.Ghi[t + 1][j];
                }
                m_mats.Ghi[k] = _mm256_load_si256(reinterpret_cast<const __m256i *>(buf));
            }

            // Константы по модулям: дорожки 0..3 — p1, дорожки 4..7 — p2.
            alignas(32) u32 p32[8], pm1[8], r32[8], odd[8];
            for (int i = 0; i < 8; ++i)
            {
                const u32 p = (i < 4) ? static_cast<u32>(p1) : static_cast<u32>(p2);
                p32[i] = p;
                pm1[i] = p - 1;
                r32[i] = static_cast<u32>(0x100000000ull / p);
            }
            for (int i = 0; i < 4; ++i)
            {
                odd[i * 2] = r32[i * 2 + 1];
                odd[i * 2 + 1] = 0;
            }
            m_mats.P = _mm256_load_si256(reinterpret_cast<const __m256i *>(p32));
            m_mats.Pm1 = _mm256_load_si256(reinterpret_cast<const __m256i *>(pm1));
            m_mats.R32 = _mm256_load_si256(reinterpret_cast<const __m256i *>(r32));
            m_mats.R32odd = _mm256_load_si256(reinterpret_cast<const __m256i *>(odd));

            alignas(32) u16 p16[16], m16[16];
            for (int i = 0; i < 16; ++i)
            {
                const u16 p = (i < 8) ? static_cast<u16>(p1) : static_cast<u16>(p2);
                p16[i] = p;
                m16[i] = static_cast<u16>(65536 / p);
            }
            m_mats.P16 = _mm256_load_si256(reinterpret_cast<const __m256i *>(p16));
            m_mats.M16 = _mm256_load_si256(reinterpret_cast<const __m256i *>(m16));

            // Одна и та же маска выбирает пару соседних 16-битных ячеек
            // (state[i], state[i+1]) в обеих 128-битных половинах.
            alignas(32) char shuf[32];
            auto mk_shuf = [&](int b)
            {
                for (int i = 0; i < 8; ++i)
                {
                    shuf[i * 4 + 0] = static_cast<char>(b + 0);
                    shuf[i * 4 + 1] = static_cast<char>(b + 1);
                    shuf[i * 4 + 2] = static_cast<char>(b + 2);
                    shuf[i * 4 + 3] = static_cast<char>(b + 3);
                }
                return _mm256_load_si256(reinterpret_cast<const __m256i *>(shuf));
            };
            m_mats.m01 = mk_shuf(0);
            m_mats.m23 = mk_shuf(4);
            m_mats.m45 = mk_shuf(8);
            m_mats.m67 = mk_shuf(12);

            alignas(32) u64 inv[4] = {0, 0, ~0ull, ~0ull};
            m_mats.hi_inv = _mm256_load_si256(reinterpret_cast<const __m256i *>(inv));
        }

        /// Приведение входа по модулю p. Результат лежит в [0, 2p), чего
        /// достаточно: финальный модуль на блоке всё равно берётся.
        __m256i reduce_p16(__m256i v) const
        {
            const __m256i q = _mm256_mulhi_epu16(v, m_mats.M16);
            return _mm256_sub_epi16(v, _mm256_mullo_epi16(q, m_mats.P16));
        }

        /// Редукция восьми 32-битных накоплений по модулю p (Барретт).
        __m256i mod_p32(__m256i x) const
        {
            const __m256i pe = _mm256_mul_epu32(x, m_mats.R32);
            const __m256i po = _mm256_mul_epu32(_mm256_srli_epi64(x, 32), m_mats.R32odd);
            const __m256i q = _mm256_blend_epi32(_mm256_srli_epi64(pe, 32),
                                                 _mm256_slli_epi64(_mm256_srli_epi64(po, 32), 32), 0xAA);
            const __m256i r = _mm256_sub_epi32(x, _mm256_mullo_epi32(q, m_mats.P));
            const __m256i mask = _mm256_cmpgt_epi32(r, m_mats.Pm1);
            return _mm256_sub_epi32(r, _mm256_and_si256(mask, m_mats.P));
        }
#else
        /// Предвычисленные векторы SSE4.2-прохода (на один генератор).
        struct mats128
        {
            __m128i A01lo, A23lo, A01hi, A23hi;
            __m128i Glo[STEPS / 2], Ghi[STEPS / 2];
            __m128i P, Pm1, R32;     // 32-битная редукция по модулю.
            __m128i P16, M16;        // 16-битная редукция входа.
            __m128i m01, m23, m45, m67; // pshufb: пара соседних ячеек.
            __m128i ff;              // 0xFFFF: инверсия входных байтов.
        };

        alignas(16) mats128 m_p1{};
        alignas(16) mats128 m_p2{};

        static void build_one(const trans &t, u32 p, mats128 &m)
        {
            alignas(16) u16 buf[8];

            auto pair = [&](int i0, int i1, const u16 (&A)[4][4])
            {
                for (int j = 0; j < 4; ++j)
                {
                    buf[j * 2] = A[i0][j];
                    buf[j * 2 + 1] = A[i1][j];
                }
                return _mm_load_si128(reinterpret_cast<const __m128i *>(buf));
            };

            m.A01lo = pair(0, 1, t.Alo);
            m.A23lo = pair(2, 3, t.Alo);
            m.A01hi = pair(0, 1, t.Ahi);
            m.A23hi = pair(2, 3, t.Ahi);

            for (int k = 0; k < STEPS / 2; ++k)
            {
                const int tt = 2 * k;
                for (int j = 0; j < 4; ++j)
                {
                    buf[j * 2] = t.Glo[tt][j];
                    buf[j * 2 + 1] = t.Glo[tt + 1][j];
                }
                m.Glo[k] = _mm_load_si128(reinterpret_cast<const __m128i *>(buf));

                for (int j = 0; j < 4; ++j)
                {
                    buf[j * 2] = t.Ghi[tt][j];
                    buf[j * 2 + 1] = t.Ghi[tt + 1][j];
                }
                m.Ghi[k] = _mm_load_si128(reinterpret_cast<const __m128i *>(buf));
            }

            m.P = _mm_set1_epi32(static_cast<int>(p));
            m.Pm1 = _mm_set1_epi32(static_cast<int>(p - 1));
            m.R32 = _mm_set1_epi32(static_cast<int>(0x100000000ull / p));
            m.P16 = _mm_set1_epi16(static_cast<short>(p));
            m.M16 = _mm_set1_epi16(static_cast<short>(65536 / p));

            // Каждая маска размножает пару соседних ячеек по всем 32-битным дорожкам.
            m.m01 = _mm_setr_epi8(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3);
            m.m23 = _mm_setr_epi8(4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7);
            m.m45 = _mm_setr_epi8(8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11, 8, 9, 10, 11);
            m.m67 = _mm_setr_epi8(12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15);
            m.ff = _mm_set1_epi16(-1);
        }

        void build(const u16x8 &K1, const u16x8 &K2)
        {
            build_one(build_trans<p1>(K1), static_cast<u32>(p1), m_p1);
            build_one(build_trans<p2>(K2), static_cast<u32>(p2), m_p2);
        }

        /// Приведение входа по модулю p (значение остаётся в [0, 2p)).
        static __m128i reduce_p16(__m128i v, const mats128 &m)
        {
            const __m128i q = _mm_mulhi_epu16(v, m.M16);
            return _mm_sub_epi16(v, _mm_mullo_epi16(q, m.P16));
        }

        /// Редукция четырёх 32-битных накоплений по модулю p (Барретт).
        static __m128i mod_p32(__m128i x, const mats128 &m)
        {
            const __m128i pe = _mm_mul_epu32(x, m.R32);
            const __m128i po = _mm_mul_epu32(_mm_srli_epi64(x, 32), m.R32);
            const __m128i q = _mm_blend_epi16(_mm_srli_epi64(pe, 32),
                                              _mm_slli_epi64(_mm_srli_epi64(po, 32), 32), 0xCC);
            const __m128i r = _mm_sub_epi32(x, _mm_mullo_epi32(q, m.P));
            const __m128i mask = _mm_cmpgt_epi32(r, m.Pm1);
            return _mm_sub_epi32(r, _mm_and_si128(mask, m.P));
        }
#endif
    };

    /**
     * @brief Слитное возведение состояния в степень сразу для двух
     * спаренных генераторов с разными модулями (p1 и p2).
     * @details
     * `gens::add_salt` возводит в степень q оба генератора с одинаковым q.
     * Раньше это были два независимых прохода, каждый со своими 128-битными
     * SIMD-операциями. Здесь состояния пары (p1) и пары (p2) лежат в одном
     * 256-битном регистре (младшая половина — p1, старшая — p2), поэтому
     * такт генератора выполняется один раз на обе пары.
     *
     * Арифметика идентична LFSR_paired_2x4::power_by: те же скалярные
     * свёртки для вычисления входа и тот же 16-битный Барретт, только с
     * разными p в половинах регистра. Результат побитно совпадает с
     * раздельным вычислением (это проверяет эталонный тест).
     */
    template <int p1, int p2>
    class LFSR_fused_exp
    {
        static_assert(p1 > 1 && p1 < 256);
        static_assert(p2 > 1 && p2 < 256);

    public:
        LFSR_fused_exp(const u16x8 &K1, const u16x8 &K2)
        {
#ifdef USE_AVX2
            m_K = _mm256_set_m128i(load128(K2), load128(K1));
            m_M = _mm256_set_m128i(_mm_set1_epi16(static_cast<short>(M(p2))),
                                   _mm_set1_epi16(static_cast<short>(M(p1))));
            m_P = _mm256_set_m128i(_mm_set1_epi16(static_cast<short>(p2)),
                                   _mm_set1_epi16(static_cast<short>(p1)));
            m_Pm1 = _mm256_set_m128i(_mm_set1_epi16(static_cast<short>(p2 - 1)),
                                     _mm_set1_epi16(static_cast<short>(p1 - 1)));

            // pshufb: размножить ячейку 3 (байты 6..7) на нижние 4 ячейки
            // и ячейку 7 (байты 14..15) — на верхние 4, в каждом lane.
            m_bcast = _mm256_setr_epi8(
                6, 7, 6, 7, 6, 7, 6, 7, 14, 15, 14, 15, 14, 15, 14, 15,
                6, 7, 6, 7, 6, 7, 6, 7, 14, 15, 14, 15, 14, 15, 14, 15);

            // Занулить ячейки 0 и 4 (байты 0..1 и 8..9) в каждом lane.
            m_iso = _mm256_setr_epi8(
                0, 0, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1,
                0, 0, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1);
#else
            m_K1 = load128(K1);
            m_K2 = load128(K2);
            m_M1 = _mm_set1_epi16(static_cast<short>(M(p1)));
            m_M2 = _mm_set1_epi16(static_cast<short>(M(p2)));
            m_P1 = _mm_set1_epi16(static_cast<short>(p1));
            m_P2 = _mm_set1_epi16(static_cast<short>(p2));
            m_Pm1_1 = _mm_set1_epi16(static_cast<short>(p1 - 1));
            m_Pm1_2 = _mm_set1_epi16(static_cast<short>(p2 - 1));
            m_bcast = _mm_setr_epi8(
                6, 7, 6, 7, 6, 7, 6, 7, 14, 15, 14, 15, 14, 15, 14, 15);
            m_iso = _mm_setr_epi8(
                0, 0, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1);
#endif
        }

        /**
         * @brief state1 := state1^q, state2 := state2^q.
         * @details q = 0 — тождество (состояния не меняются).
         */
        void run(u16x8 &state1, u16x8 &state2, u64 q) const
        {
            if (q == 0)
                return;

            u16x8 acc1{1, 0, 0, 0, 1, 0, 0, 0};
            u16x8 acc2 = acc1;
            u16x8 src1 = state1;
            u16x8 src2 = state2;

            for (;;)
            {
                if ((q & 1) != 0)
                    mul(acc1, src1, acc2, src2);
                q >>= 1;
                if (q == 0)
                    break; // последний квадрат уже не нужен
                mul(src1, src1, src2, src2);
            }

            state1 = acc1;
            state2 = acc2;
        }

    private:
        static constexpr u16 M(int p) { return static_cast<u16>((1u << 16) / p); }

        static __m128i load128(const u16x8 &a)
        {
            // Локальные копии состояний (acc/src) не имеют alignas(16),
            // поэтому используем невыровненную загрузку.
            return _mm_loadu_si128(reinterpret_cast<const __m128i *>(a.data()));
        }

        /// acc := acc * b для обеих пар; b может совпадать с acc (square).
        void mul(u16x8 &a1, const u16x8 &b1, u16x8 &a2, const u16x8 &b2) const
        {
            u16x8 ac1 = a1, ac2 = a2;
            const u16x8 &br1 = (&b1 == &a1) ? ac1 : b1;
            const u16x8 &br2 = (&b2 == &a2) ? ac2 : b2;

            // Как и в mult_by, накопитель начинает с нуля: произведение
            // собирается семью тактами next() из скалярных свёрток.
#ifdef USE_AVX2
            __m256i st = _mm256_setzero_si256();
#else
            __m128i st1 = _mm_setzero_si128(), st2 = _mm_setzero_si128();
#endif

            for (int power = 6; power >= 0; --power)
            {
                u32 v1 = 0, v2 = 0, v3 = 0, v4 = 0;
                const int first = (power > 3) ? power - 3 : 0;
                const int last = (power < 3) ? power : 3;
                for (int i = first; i <= last; ++i)
                {
                    v1 += static_cast<u32>(ac1[i]) * br1[power - i];
                    v2 += static_cast<u32>(ac1[i + 4]) * br1[power - i + 4];
                    v3 += static_cast<u32>(ac2[i]) * br2[power - i];
                    v4 += static_cast<u32>(ac2[i + 4]) * br2[power - i + 4];
                }

#ifdef USE_AVX2
                st = next(st, static_cast<u16>(v1 % p1), static_cast<u16>(v2 % p1),
                              static_cast<u16>(v3 % p2), static_cast<u16>(v4 % p2));
#else
                st1 = next(st1, m_K1, m_M1, m_P1, m_Pm1_1,
                           static_cast<u16>(v1 % p1), static_cast<u16>(v2 % p1));
                st2 = next(st2, m_K2, m_M2, m_P2, m_Pm1_2,
                           static_cast<u16>(v3 % p2), static_cast<u16>(v4 % p2));
#endif
            }

#ifdef USE_AVX2
            _mm_storeu_si128(reinterpret_cast<__m128i *>(a1.data()),
                             _mm256_castsi256_si128(st));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(a2.data()),
                             _mm256_extracti128_si256(st, 1));
#else
            _mm_storeu_si128(reinterpret_cast<__m128i *>(a1.data()), st1);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(a2.data()), st2);
#endif
        }

#ifdef USE_AVX2
        /// Один такт обеих пар: сдвиг, обратная связь K, входы, Барретт.
        __m256i next(__m256i st, u16 v1, u16 v2, u16 v3, u16 v4) const
        {
            const __m256i shifted = _mm256_and_si256(_mm256_slli_si256(st, 2), m_iso);
            const __m256i vv = _mm256_shuffle_epi8(st, m_bcast);
            const __m256i prod = _mm256_mullo_epi16(vv, m_K);
            const __m256i inp = _mm256_setr_epi16(
                static_cast<short>(v1), 0, 0, 0, static_cast<short>(v2), 0, 0, 0,
                static_cast<short>(v3), 0, 0, 0, static_cast<short>(v4), 0, 0, 0);

            __m256i res = _mm256_add_epi16(_mm256_add_epi16(shifted, prod), inp);
            const __m256i q = _mm256_mulhi_epu16(res, m_M);
            res = _mm256_sub_epi16(res, _mm256_mullo_epi16(q, m_P));
            const __m256i mask = _mm256_cmpgt_epi16(res, m_Pm1);
            return _mm256_sub_epi16(res, _mm256_and_si256(mask, m_P));
        }
#else
        /// Один такт одной пары (для SSE4.2-сборки).
        __m128i next(__m128i st, __m128i K, __m128i Mv, __m128i Pv, __m128i Pm1,
                     u16 v1, u16 v2) const
        {
            const __m128i shifted = _mm_and_si128(_mm_slli_si128(st, 2), m_iso);
            const __m128i vv = _mm_shuffle_epi8(st, m_bcast);
            const __m128i prod = _mm_mullo_epi16(vv, K);
            const __m128i inp = _mm_setr_epi16(
                static_cast<short>(v1), 0, 0, 0, static_cast<short>(v2), 0, 0, 0);

            __m128i res = _mm_add_epi16(_mm_add_epi16(shifted, prod), inp);
            const __m128i q = _mm_mulhi_epu16(res, Mv);
            res = _mm_sub_epi16(res, _mm_mullo_epi16(q, Pv));
            const __m128i mask = _mm_cmpgt_epi16(res, Pm1);
            return _mm_sub_epi16(res, _mm_and_si128(mask, Pv));
        }
#endif

#ifdef USE_AVX2
        __m256i m_K{}, m_M{}, m_P{}, m_Pm1{}, m_bcast{}, m_iso{};
#else
        __m128i m_K1{}, m_K2{}, m_M1{}, m_M2{}, m_P1{}, m_P2{};
        __m128i m_Pm1_1{}, m_Pm1_2{}, m_bcast{}, m_iso{};
#endif
    };

}