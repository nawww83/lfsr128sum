#pragma once

#include "lfsr.h"
#include <algorithm>
#include <bit>
#include <cstring>
#include <numeric>
#include <span>
#include <stdexcept>
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
        int q;  // Количество тактов, применяемое к спаренному генератору.
                // q = 0 означает «без сдвига»: возведение в степень не применяется
                // (иначе x^0 = 1 стирало бы состояние).
        u16 s0; // Константный вход первого генератора.
        u16 s1; // -//- второго генератора.
    };

    // Периоды четырёх регистров. Показатель соли q обязан быть взаимно прост
    // с каждым из них: x -> x^q — гомоморфизм с ядром gcd(q, T), поэтому при
    // gcd > 1 образ состояния сжимается ровно в gcd(q, T) раз (и вместе с ним
    // падает эффективное пространство состояний, а значит и birthday-граница).
    static constexpr u64 T_251_4 = 251ull * 251 * 251 * 251 - 1;
    static constexpr u64 T_251_3 = 251ull * 251 * 251 - 1;
    static constexpr u64 T_241_4 = 241ull * 241 * 241 * 241 - 1;
    static constexpr u64 T_241_3 = 241ull * 241 * 241 - 1;

    /// Инъективен ли показатель соли на всех четырёх регистрах.
    constexpr bool salt_is_injective(int q)
    {
        return (q > 0) &&
               std::gcd(static_cast<u64>(q), T_251_4) == 1 &&
               std::gcd(static_cast<u64>(q), T_251_3) == 1 &&
               std::gcd(static_cast<u64>(q), T_241_4) == 1 &&
               std::gcd(static_cast<u64>(q), T_241_3) == 1;
    }

    static constexpr salt S0{61, 2, 3};   // Доменная соль старта: считается один раз на процесс.
    static constexpr salt S1{13, 4, 7};   // Первая соль финализации.
    static constexpr salt S2{37, 2, 13};  // Вторая соль финализации.

    static_assert(S0.q >= 4); // Минимально достаточное насыщение, m = 4
    static_assert(S1.q >= 4);
    static_assert(S2.q >= 4);

    static_assert(salt_is_injective(S0.q));
    static_assert(salt_is_injective(S1.q));
    static_assert(salt_is_injective(S2.q));

    struct gens
    {
        LFSR251x4 g_251x4; // Два спаренных генератора с периодами T1 = p1^4 - 1, T2 = p1^3 - 1: НОД(T1, T2) = p1 - 1.
        LFSR241x4 g_241x4; // Два спаренных генератора с периодами T1 = p2^4 - 1, T2 = p2^3 - 1: НОД(T1, T2) = p2 - 1.
        // Простые числа p1, p2 подобраны так, чтобы минимизировать НОД, но и стараться не выходить за пределы байта.

    public:
        /**
         * @brief Общий слитный блоковый процессор (AVX2 или SSE4.2).
         * @details Матрицы зависят только от K1/K2 и модулей, поэтому строятся
         * один раз на процесс; инициализация статики потокобезопасна.
         * Реализация выбирается на этапе компиляции, результаты совпадают.
         */
        static const lfsr8::LFSR_fused_2x4<251, 241> &fused_engine()
        {
            static const lfsr8::LFSR_fused_2x4<251, 241> engine{K1, K2};
            return engine;
        }

        /**
         * @brief Готовые матрицы первого генератора (только для чтения).
         * @details Матрицы зависят только от K1/K2, поэтому считаются один раз
         * на процесс и копируются в каждый gens. Пересчёт в конструкторе стоил
         * ~1 мкс и был основной стоимостью хэширования мелких входов.
         */
        static const LFSR251x4 &prototype_251()
        {
            static const LFSR251x4 engine{K1};
            return engine;
        }

        /// Готовые матрицы второго генератора (только для чтения).
        static const LFSR241x4 &prototype_241()
        {
            static const LFSR241x4 engine{K2};
            return engine;
        }

        /// Слитное возведение в степень для обеих пар (251 и 241).
        static const lfsr8::LFSR_fused_exp<251, 241> &fused_exp()
        {
            static const lfsr8::LFSR_fused_exp<251, 241> engine{K1, K2};
            return engine;
        }

        /**
         * @brief Конструктор по умолчанию.
         * @details Копирует предвычисленные матрицы из статических прототипов
         * (единицы наносекунд) вместо повторного пересчёта.
         */
        gens() : g_251x4(prototype_251()),
                 g_241x4(prototype_241()) {}
        
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
         * При q = 0 возведение в степень не применяется (состояние не стирается):
         * так соль с нулевым сдвигом остаётся безобидной.
         */
        inline void add_salt(const salt& S)
        {
            g_251x4.next_simd(S.s0);
            g_241x4.next_simd(S.s1);

            // Возведение в степень одинаково для обеих пар, поэтому делаем
            // его слитно в одном 256-битном регистре. Если коэффициенты
            // подменили через set_K(), слитный путь неприменим.
            if (g_251x4.get_generator_coeffs() == K1 &&
                g_241x4.get_generator_coeffs() == K2)
            {
                fused_exp().run(g_251x4.state_ref(), g_241x4.state_ref(), S.q);
            }
            else
            {
                g_251x4.power_by(S.q);
                g_241x4.power_by(S.q);
            }
        }

        /**
         * @brief Состояние после доменной соли S0, посчитанное один раз.
         * @details Стартовая соль применяется всегда к одной и той же точке
         * (единичному состоянию), поэтому возведение в степень незачем
         * повторять на каждый хэш: прототип считается при первом обращении
         * и дальше только копируется. Это снимает фиксированную стоимость
         * с каждого вызова, включая самые мелкие входы.
         */
        static const gens &start()
        {
            static const gens proto = []
            {
                gens g;
                g.add_salt(S0);
                return g;
            }();
            return proto;
        }
    };

    /**
     * @brief Довести остаток сообщения (< размера слитного блока) до состояния.
     * @details Полные 16-байтные блоки идут через матричный SIMD-путь, хвост
     * (< 16 байт) — побайтово с двух концов. Используется hasher::finish()
     * и тестами; отдельного поблочного API у хэшера нет.
     */
    inline void process_remainder(gens &g, std::span<const std::byte> input)
    {
        const size_t n = input.size();
        if (n == 0)
            return;

        size_t i = 0;
        for (; i + 16 <= n; i += 16)
        {
            const __m128i data = _mm_loadu_si128(
                reinterpret_cast<const __m128i *>(input.data() + i));
            g.g_251x4.next_simd_block(data);
            g.g_241x4.next_simd_block(_mm_xor_si128(data, _mm_set1_epi8(0xFF)));
        }

        if (i < n)
        {
            size_t j = n;
            while (i < j)
            {
                const uint16_t l = static_cast<uint8_t>(input[i]);
                const uint16_t r = (i + 1 < j) ? static_cast<uint8_t>(input[j - 1]) : 0;
                g.g_251x4.next_simd(l, r);
                g.g_241x4.next_simd(static_cast<uint16_t>(~r), static_cast<uint16_t>(~l));
                ++i;
                if (j > i)
                    --j;
            }
        }
    }

    // ======================================================================
    // Потоковый хэшер с нелинейным слоем.
    //
    // Зачем нелинейный слой. Инъекция входа в LFSR ЛИНЕЙНА над GF(p): в
    // состояние попадает только значение слова по модулю p, умноженное на
    // постоянную матрицу перехода. Поэтому два разных 16-битных слова с
    // одинаковым остатком дают одинаковое состояние, а значит и одинаковый
    // хэш: для p1=251 и p2=241 достаточно сдвига Δ = lcm(251,241) = 60491,
    // чтобы обе ветки «не заметили» правку. Это превращает любую линейную
    // схему (CRC, LFSR) в конструктор бесплатных коллизий.
    //
    // Нелинейный ARX-фолд добавляет к состоянию 128-битное слагаемое,
    // которое зависит от данных через умножения и переносы. Теперь правка
    // данных обязана одновременно обнулить и линейный вклад, и фолд, что
    // сводит задачу к поиску коллизии 128-битного значения.
    //
    // Раскладка фолда: 4 цепочки по 32 байта (8 x u32). Каждая цепочка
    // обновляется раз в 128 байт, поэтому латентность умножения не упирается
    // в конвейер; раз в 512 байт цепочки связываются между собой, чтобы
    // разностный след не мог остаться в одной полосе. Результат не зависит
    // ни от размера кусков, подаваемых в update(), ни от наличия AVX2:
    // SSE-путь считает ту же математику в двух 128-битных половинах.
    // ======================================================================

    /// Вектор из 8 x u32: 256 бит в AVX2-сборке, две 128-битные половины иначе.
    struct lane8
    {
#ifdef USE_AVX2
        __m256i v;
#else
        __m128i lo;
        __m128i hi;
#endif
    };

    static constexpr u32 FOLD_PRIME = 0x9E3779B1u;

    inline void fold_load(lane8 &d, const uint8_t *p)
    {
#ifdef USE_AVX2
        d.v = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
#else
        d.lo = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
        d.hi = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p + 16));
#endif
    }

    inline void fold_add(lane8 &a, const lane8 &b)
    {
#ifdef USE_AVX2
        a.v = _mm256_add_epi32(a.v, b.v);
#else
        a.lo = _mm_add_epi32(a.lo, b.lo);
        a.hi = _mm_add_epi32(a.hi, b.hi);
#endif
    }

    inline void fold_mul(lane8 &a)
    {
#ifdef USE_AVX2
        a.v = _mm256_mullo_epi32(a.v, _mm256_set1_epi32(FOLD_PRIME));
#else
        const __m128i c = _mm_set1_epi32(static_cast<int>(FOLD_PRIME));
        a.lo = _mm_mullo_epi32(a.lo, c);
        a.hi = _mm_mullo_epi32(a.hi, c);
#endif
    }

    inline void fold_rotl13(lane8 &a)
    {
#ifdef USE_AVX2
        a.v = _mm256_or_si256(_mm256_slli_epi32(a.v, 13), _mm256_srli_epi32(a.v, 19));
#else
        a.lo = _mm_or_si128(_mm_slli_epi32(a.lo, 13), _mm_srli_epi32(a.lo, 19));
        a.hi = _mm_or_si128(_mm_slli_epi32(a.hi, 13), _mm_srli_epi32(a.hi, 19));
#endif
    }

    /// Поменять местами 32-битные половины каждого u32 (одна инструкция).
    inline void fold_shuf(lane8 &a)
    {
#ifdef USE_AVX2
        a.v = _mm256_shuffle_epi32(a.v, _MM_SHUFFLE(2, 3, 0, 1));
#else
        a.lo = _mm_shuffle_epi32(a.lo, _MM_SHUFFLE(2, 3, 0, 1));
        a.hi = _mm_shuffle_epi32(a.hi, _MM_SHUFFLE(2, 3, 0, 1));
#endif
    }

    /**
     * @brief Один раунд ARX над 8 x u32 (модифицированный раунд ChaCha).
     * @details Общий для нелинейного фолда и финализации: одна и та же
     * перемешивающая функция с разным числом раундов.
     */
    inline void arx_round8(u32 w[8])
    {
        for (int i = 0; i < 4; ++i)
        {
            u32 &x = w[i], &y = w[i + 4];
            u32 &z = w[(i + 1) & 3], &t = w[4 + ((i + 1) & 3)];
            x += y; t ^= x; t = std::rotl(t, 16);
            z += t; y ^= z; y = std::rotl(y, 12);
            x += y; t ^= x; t = std::rotl(t, 8);
            z += t; y ^= z; y = std::rotl(y, 7);
        }
    }

    inline void fold_store(const lane8 &a, u32 out[8])
    {
#ifdef USE_AVX2
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(out), a.v);
#else
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out), a.lo);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out + 4), a.hi);
#endif
    }

    /**
     * @brief Нелинейный ARX-фолд сообщения.
     * @details Потоковый: состояние собирается по кускам любого размера и
     * зависит только от порядка байт, но не от границ вызовов update().
     */
    class fold_acc
    {
    public:
        /**
         * @brief Сброс аккумулятора.
         * @param seed 8 слов начального состояния; nullptr — константный сид.
         * @details В MAC-режиме сид зависит от ключа, поэтому нелинейный
         * слой отличается с первого блока, а не только в финализации.
         */
        void reset(const u32 *seed = nullptr)
        {
            static constexpr u32 default_seed[8] = {0x9E3779B1u, 0x85EBCA77u, 0xC2B2AE3Du, 0x27D4EB2Fu,
                                                    0x165667B1u, 0xD3A2646Cu, 0xFD7046C5u, 0xB55A4F09u};
            if (!seed)
                seed = default_seed;
            for (int i = 0; i < 4; ++i)
            {
                u32 w[8];
                for (int j = 0; j < 8; ++j)
                    w[j] = seed[j] + static_cast<u32>(i) * 0x01010101u;
                fold_load(m_lane[i], reinterpret_cast<const uint8_t *>(w));
            }
            m_pos = 0;
            m_tail_len = 0;
        }

        /// Свернуть очередные байты. Внутренняя буферизация хвоста < 32 байт.
        void update(const uint8_t *p, size_t n)
        {
            if (m_tail_len)
            {
                const size_t take = (std::min)(n, 32 - m_tail_len);
                std::memcpy(m_tail + m_tail_len, p, take);
                m_tail_len += take;
                p += take;
                n -= take;
                if (m_tail_len == 32)
                {
                    fold32(m_tail);
                    m_tail_len = 0;
                }
            }
            while (n >= 32)
            {
                fold32(p);
                p += 32;
                n -= 32;
            }
            if (n)
            {
                std::memcpy(m_tail, p, n);
                m_tail_len = n;
            }
        }

        /// Довернуть нулями незакрытый хвост (детерминированное завершение).
        void flush()
        {
            if (m_tail_len)
            {
                std::memset(m_tail + m_tail_len, 0, 32 - m_tail_len);
                fold32(m_tail);
                m_tail_len = 0;
            }
        }

        /**
         * @brief Финализация: 128 бит из фолда с привязкой полной длины.
         * @param total_len Полная длина сообщения в байтах.
         * @param out 16 байт нелинейного слагаемого.
         */
        void finalize(u64 total_len, uint8_t out[16]) const
        {
            u32 w[8];
            fold_store(m_lane[0], w);
            for (int i = 1; i < 4; ++i)
            {
                u32 t[8];
                fold_store(m_lane[i], t);
                for (int j = 0; j < 8; ++j)
                    w[j] ^= t[j];
            }

            // Привязка ПОЛНОЙ длины (не только чётности и не только % 64 КиБ).
            const u32 c0 = w[0], c5 = w[5];
            w[0] ^= static_cast<u32>(total_len);
            w[1] ^= static_cast<u32>(total_len >> 32);
            w[2] += static_cast<u32>(total_len * 0x9E3779B1ull);
            w[3] += static_cast<u32>((total_len * 0xC2B2AE3Dull) >> 17);
            w[4] += static_cast<u32>(total_len * 0x27D4EB2Full);
            w[5] ^= static_cast<u32>((total_len * 0x85EBCA77ull) >> 29);
            w[6] += c0 + w[5];
            w[7] ^= w[1] + c5;

            // Три раунда в духе ChaCha: add / xor / rotl (нелинейно и над
            // GF(2), и над GF(p) из-за переносов).
            for (int r = 0; r < 3; ++r)
                arx_round8(w);
            for (int j = 0; j < 8; ++j)
            {
                out[2 * j] = static_cast<uint8_t>(w[j] ^ (w[j] >> 16));
                out[2 * j + 1] = static_cast<uint8_t>((w[j] << 7) ^ (w[j] >> 9));
            }
        }

        /// Аккумулятор полосы (только для чтения; используется тестами).
        const lane8 &lane(int i) const
        {
            return m_lane[i];
        }

    private:
        void fold32(const uint8_t *p)
        {
            const int lane = static_cast<int>((m_pos / 32) & 3);
            lane8 v;
            fold_load(v, p);
            fold_add(m_lane[lane], v);
            fold_mul(m_lane[lane]);
            fold_rotl13(m_lane[lane]);
            m_pos += 32;

            // Раз в 512 байт полосы связываются между собой: четыре
            // «перестановки половин слова + сложение» по кругу. Сложение
            // (а не XOR) выбрано специально: XOR по замкнутому циклу из
            // четырёх шагов с одной и той же инволюцией гасил бы разность
            // (δ ⊕ δ = 0), а переносы плюс разный порядок шагов — нет.
            // Это дешёвая (≈2% потока) линейная кросс-связь, не нелинейный
            // барьер: нелинейность даёт умножение фолда, а барьер против
            // подгонки — финализация и линейная часть.
            if ((m_pos & 511) == 0)
            {
                lane8 t;
                t = m_lane[0]; fold_shuf(t); fold_add(m_lane[1], t);
                t = m_lane[1]; fold_shuf(t); fold_add(m_lane[2], t);
                t = m_lane[2]; fold_shuf(t); fold_add(m_lane[3], t);
                t = m_lane[3]; fold_shuf(t); fold_add(m_lane[0], t);
            }
        }

        lane8 m_lane[4];
        u64 m_pos = 0;
        uint8_t m_tail[32] = {};
        size_t m_tail_len = 0;
    };

    // ======================================================================
    // Keyed-режим (MAC).
    //
    // Ключ входит в три места сразу:
    //   1) стартовое состояние — ключ прогоняется тем же линейным движком,
    //      что и сообщение, затем применяется доменная соль S0;
    //   2) сид нелинейного фолда — нелинейный слой отличается с первого байта;
    //   3) слова финализации — XOR до feed-forward, поэтому разность ключа
    //      проходит через пять раундов ARX, а не складывается с выходом.
    //
    // Почему не HMAC-обёртка. Итеративной цепочки «выход блока -> вход
    // следующего» здесь нет (сообщение входит в состояние целиком), поэтому
    // length extension в MD-смысле неприменим, и один проход дешевле двух
    // в HMAC. Обратная сторона: это домашняя конструкция, а не проверенная
    // связка, и на стойкость MAC как PRF здесь никто не заявляет.
    // ======================================================================

    /// Максимальная длина ключа. Длиннее — это уже не ключ, а пароль:
    /// его надо прогнать через KDF (argon2/scrypt) и взять 32 байта.
    inline constexpr size_t KEY_MAX = 64;

    /// Стартовое состояние: 16 ячеек обоих генераторов.
    struct start_pt
    {
        STATE s251;
        STATE s241;
    };

    /**
     * @brief Предвычисленный ключевой контекст.
     * @details Считается один раз на ключ и переиспользуется: на каждое
     * сообщение остаётся только копия состояния, как и в unkeyed-режиме.
     * Поэтому keyed-режим не стоит на сообщение ничего.
     */
    struct key_ctx
    {
        start_pt start;   // стартовое состояние (ключ + доменная соль)
        u32 fold_seed[8]; // ключевой сид нелинейного фолда
        u32 kw[8];        // ключевые слова финализации
    };

    /// Развести 8 слов: 6 раундов ARX + feed-forward (как в финализации).
    inline void key_mix(u32 w[8])
    {
        u32 ff[8];
        std::memcpy(ff, w, sizeof(ff));
        for (int r = 0; r < 6; ++r)
            arx_round8(w);
        for (int i = 0; i < 8; ++i)
            w[i] ^= ff[i];
    }

    /**
     * @brief Построить ключевой контекст.
     * @param key Ключ длиной 1..KEY_MAX байт.
     * @throws std::invalid_argument при пустом или слишком длинном ключе.
     * @details Ключ дополняется нулями до 64 байт (как в HMAC), а его длина
     * входит в доменную метку: ключи "a" и "a\0" дают разные теги.
     */
    inline key_ctx make_key(std::span<const std::byte> key)
    {
        if (key.empty() || key.size() > KEY_MAX)
            throw std::invalid_argument("lfsr128 MAC: длина ключа должна быть 1..64 байта");

        uint8_t kb[KEY_MAX] = {};
        std::memcpy(kb, key.data(), key.size());

        uint8_t tag[16] = {};
        static constexpr char DOM[] = "lfsr128-mac-v4"; // 15 символов + длина ключа
        std::memcpy(tag, DOM, sizeof(DOM) - 1);
        tag[15] = static_cast<uint8_t>(key.size());

        key_ctx kc;

        // Стартовое состояние: доменная метка и ключ проходят через тот же
        // линейный движок, что и сообщение, затем применяется соль S0.
        gens g;
        g.reset();
        process_remainder(g, std::as_bytes(std::span(tag, sizeof(tag))));
        process_remainder(g, std::as_bytes(std::span(kb, sizeof(kb))));
        g.add_salt(S0);
        kc.start = {g.g_251x4.get_state(), g.g_241x4.get_state()};

        // Ключевые слова: 16 слов ключевого блока сворачиваются в 8 двумя
        // разными свёртками и разводятся ARX.
        u32 w[16];
        std::memcpy(w, kb, sizeof(w));
        for (int i = 0; i < 8; ++i)
        {
            kc.kw[i] = w[i] ^ w[i + 8];
            kc.fold_seed[i] = w[i] + w[i + 8];
        }
        key_mix(kc.kw);
        key_mix(kc.fold_seed);
        return kc;
    }

    /**
     * @brief Потоковый 128-битный хэшер v4.
     * @details Один проход по данным: линейный LFSR-движок + нелинейный фолд,
     * результат фолда (с полной длиной) инжектируется в состояние перед
     * финализацией. Данные можно подавать кусками любого размера.
     *
     * Скорость на i7-12700: ~7.9 ГБ/с в файловом режиме. Соль применяется
     * один раз при старте (готовый прототип) и дважды при финализации, а не
     * на каждый блок; нелинейный фолд стоит около 17% на одном большом
     * вызове.
     */
    class hasher
    {
    public:
        hasher() { reset(); }

        /// Keyed-режим: тот же хэш, но с секретным ключом (MAC).
        explicit hasher(const key_ctx &k) { set_key(k); }

        /**
         * @brief Включить keyed-режим (MAC).
         * @details Вызывать до первого update(): контекст задаёт стартовое
         * состояние, поэтому посреди сообщения его менять нельзя.
         */
        void set_key(const key_ctx &k)
        {
            std::memcpy(m_kw, k.kw, sizeof(m_kw));
            std::memcpy(m_fold_seed, k.fold_seed, sizeof(m_fold_seed));
            m_start = k.start;
            m_keyed = true;
            reset();
        }

        /// Работает ли хэшер в режиме MAC.
        bool keyed() const { return m_keyed; }

        /// Начать новое сообщение. Ставит стартовое состояние (или S0).
        void reset()
        {
            // S0 применяется к единичному состоянию, то есть всегда к одной
            // и той же точке, — берём готовый результат вместо возведения
            // в степень на каждый вызов. В keyed-режиме стартовое состояние
            // тоже предвычислено (один раз на ключ).
            const start_pt &s = m_keyed ? m_start : unkeyed_start();
            m_g.g_251x4.set_state(s.s251);
            m_g.g_241x4.set_state(s.s241);
            m_fold.reset(m_keyed ? m_fold_seed : nullptr);
            m_len = 0;
            m_staged = 0;
        }

        /// Добавить очередной кусок данных (состояние сохраняется между вызовами).
        void update(std::span<const std::byte> data)
        {
            const uint8_t *p = reinterpret_cast<const uint8_t *>(data.data());
            size_t n = data.size();
            m_len += n;

            if (m_staged)
            {
                const size_t take = (std::min)(n, STAGE - m_staged);
                std::memcpy(m_stage + m_staged, p, take);
                m_staged += take;
                p += take;
                n -= take;
                if (m_staged == STAGE)
                {
                    m_fold.update(m_stage, STAGE);
                    run_block(m_stage, STAGE);
                    m_staged = 0;
                }
            }
            while (n >= STAGE)
            {
                m_fold.update(p, STAGE);
                run_block(p, STAGE);
                p += STAGE;
                n -= STAGE;
            }
            if (n)
            {
                std::memcpy(m_stage, p, n);
                m_staged = n;
            }
        }

        /// Завершить сообщение и вернуть 128-битный хэш.
        u128 finish()
        {
            // 1. Остаток линейного движка: полные 16-байтные блоки и хвост.
            process_remainder(m_g, std::as_bytes(std::span(m_stage, m_staged)));

            // 2. Нелинейный фолд + полная длина, инъекция в состояние.
            m_fold.update(m_stage, m_staged); // хвост, не попавший в update()
            m_fold.flush();
            uint8_t d[16];
            m_fold.finalize(m_len, d);
            inject(d);

            // 3. Соли-перестановки состояния и одна 128-битная финализация.
            m_g.add_salt(S1);
            m_g.add_salt(S2);
            return final_128(d);
        }

    private:
        /// Стартовое состояние unkeyed-режима (единичное состояние + S0).
        static const start_pt &unkeyed_start()
        {
            static const start_pt p = []
            {
                const gens &g = gens::start();
                return start_pt{g.g_251x4.get_state(), g.g_241x4.get_state()};
            }();
            return p;
        }

        /**
         * @brief Единая 128-битная финализация.
         * @details 16 ячеек состояния (по 8 бит, всего 128 бит) упаковываются
         * в четыре u32, к ним добавляются все 128 бит значения фолда — итого
         * 8 x u32. Пять раундов ARX, feed-forward и сжатие в 128 бит: одна
         * функция перемешивания на 256-битном входе, а не два независимых
         * 64-битных результата на одном состоянии.
         */
        u128 final_128(const uint8_t d[16]) const
        {
            auto s1 = m_g.g_251x4.get_state();
            auto s2 = m_g.g_241x4.get_state();

            u32 w[8];
            w[0] = static_cast<u32>(s1[0]) | (static_cast<u32>(s1[1]) << 8) |
                   (static_cast<u32>(s1[4]) << 16) | (static_cast<u32>(s1[5]) << 24);
            w[1] = static_cast<u32>(s1[2]) | (static_cast<u32>(s1[3]) << 8) |
                   (static_cast<u32>(s1[6]) << 16) | (static_cast<u32>(s1[7]) << 24);
            w[2] = static_cast<u32>(s2[0]) | (static_cast<u32>(s2[1]) << 8) |
                   (static_cast<u32>(s2[4]) << 16) | (static_cast<u32>(s2[5]) << 24);
            w[3] = static_cast<u32>(s2[2]) | (static_cast<u32>(s2[3]) << 8) |
                   (static_cast<u32>(s2[6]) << 16) | (static_cast<u32>(s2[7]) << 24);

            // Значение фолда целиком, включая ту часть, что уже ушла
            // в состояние через inject(): так все 128 бит фолда связаны
            // с выходом напрямую.
            std::memcpy(w + 4, d, 16);

            // Ключевые слова входят до feed-forward: разность ключа проходит
            // через пять раундов ARX, а не складывается с готовым выходом.
            if (m_keyed)
                for (int i = 0; i < 8; ++i)
                    w[i] ^= m_kw[i];

            u32 ff[8];
            std::memcpy(ff, w, sizeof(w));
            for (int r = 0; r < 5; ++r)
                arx_round8(w);
            for (int i = 0; i < 8; ++i)
                w[i] ^= ff[i];

            for (int i = 0; i < 4; ++i)
                w[i] ^= w[i + 4];
            return {(static_cast<u64>(w[0]) << 32) | w[1],
                    (static_cast<u64>(w[2]) << 32) | w[3]};
        }

        /// Размер порции, обрабатываемой за один вызов движка.
        /// Один и тот же в обеих сборках, чтобы результат не зависел от
        /// набора инструкций и чтобы накладные расходы update() не платились
        /// на каждые 16 байт.
        static constexpr size_t STAGE = 512;

        void run_block(const uint8_t *p, size_t n)
        {
            m_g.fused_engine().process_blocks(
                reinterpret_cast<const std::byte *>(p), n / 512,
                m_g.g_251x4.state_ref(), m_g.g_241x4.state_ref());
        }

        /// Нелинейное слагаемое в оба спаренных регистра (mod p).
        void inject(const uint8_t d[16])
        {
            auto &s1 = m_g.g_251x4.state_ref();
            auto &s2 = m_g.g_241x4.state_ref();
            for (int j = 0; j < 4; ++j)
            {
                const u32 lo = d[2 * j] + 256u * d[2 * j + 1];
                const u32 hi = d[8 + 2 * j] + 256u * d[9 + 2 * j];
                s1[j] = static_cast<u16>((s1[j] + lo) % 251);
                s1[4 + j] = static_cast<u16>((s1[4 + j] + hi) % 251);
                s2[j] = static_cast<u16>((s2[j] + hi) % 241);
                s2[4 + j] = static_cast<u16>((s2[4 + j] + lo) % 241);
            }
        }

        gens m_g;
        fold_acc m_fold;
        u64 m_len = 0;
        alignas(32) uint8_t m_stage[STAGE] = {};
        size_t m_staged = 0;

        // Keyed-режим (MAC): всё предвычислено в key_ctx, здесь только копия.
        bool m_keyed = false;
        start_pt m_start{};
        u32 m_fold_seed[8] = {};
        u32 m_kw[8] = {};
    };

    /// Удобная обёртка: хэш одного буфера целиком.
    inline u128 hash128(std::span<const std::byte> data)
    {
        hasher h;
        h.update(data);
        return h.finish();
    }

    /// Удобная обёртка: 128-битный MAC одного буфера целиком.
    inline u128 mac128(const key_ctx &k, std::span<const std::byte> data)
    {
        hasher h(k);
        h.update(data);
        return h.finish();
    }
}