#pragma once

#include <cstdint>
#include <cassert>
#include <array>
#include <type_traits>

#if defined(__x86_64__) || defined(_M_X64)
#define USE_SSE
#endif

#ifdef USE_SSE
#include <immintrin.h>
#include <smmintrin.h>
#endif

namespace lfsr8 {
using u64 = uint64_t;
using u32 = uint32_t;
using u16 = uint16_t;
using u16x8 = std::array<u16, 8>;
using u32x4 = std::array<u32, 4>;

template <int m>
class MType {
public:
    typedef typename std::conditional<(m <= 4), u32, u16>::type SAMPLE;
    typedef typename std::conditional<(m <= 4), u32x4, u16x8>::type STATE;
};

template <int m, typename T = typename MType<m>::STATE>
static auto is_zero(T st) {
    bool zf = true;
    for (int i=0; i<m; ++i) {
        zf &= (st[i] == 0);
    }
    return zf;
}

template <int p, int m>
class LFSR {
    using STATE = typename MType<m>::STATE;
    using SAMPLE = typename MType<m>::SAMPLE;
public:
    constexpr LFSR(STATE K): m_K(K) {
        static_assert(m <= 8);
        static_assert(m > 0);
        if constexpr (m > 4) {
            static_assert(p < 256);
        } else {
            static_assert(p < 256*256);
        }
        static_assert(p > 1);
        m_calc_inv_K0();
    };
    void set_state(STATE st) {
        m_state = st;
    }
    void set_unit_state() {
        m_state = {1};
    }
    void set_K(STATE K) {
        m_K = K;
        m_calc_inv_K0();
    }
    void next(SAMPLE input=0) {
#ifdef USE_SSE
        if constexpr (m > 4) {
            __m128i a = _mm_set1_epi16(m_state[m-1]);
            __m128i b = _mm_load_si128((const __m128i*)&m_K[0]);
            __m128i tmp = _mm_mullo_epi16(a, b);
            __m128i c = _mm_load_si128((const __m128i*)&m_state[0]);
            c = _mm_slli_si128(c, 2);
            __m128i mask = _mm_slli_si128(_mm_set1_epi16(-1), 2);
            __m128i inp = _mm_andnot_si128( mask, _mm_set1_epi16(input) );
            c = _mm_add_epi16(c, tmp);
            c = _mm_add_epi16(c, inp);
            _mm_store_si128((__m128i*)&m_state[0], c);
            for (int i=0; i<m; ++i) {
                m_state[i] %= p;
            }
        } else {
            __m128i a = _mm_set1_epi32(m_state[m-1]);
            __m128i b = _mm_load_si128((const __m128i*)&m_K[0]);
            __m128i tmp = _mm_mullo_epi32(a, b);
            __m128i c = _mm_load_si128((const __m128i*)&m_state[0]);
            c = _mm_slli_si128(c, 4);
            __m128i mask = _mm_slli_si128(_mm_set1_epi32(-1), 4);
            __m128i inp = _mm_andnot_si128( mask, _mm_set1_epi32(input) );
            c = _mm_add_epi32(c, tmp);
            c = _mm_add_epi32(c, inp);
            _mm_store_si128((__m128i*)&m_state[0], c);
            for (int i=0; i<m; ++i) {
                m_state[i] %= p;
            }
        }
#else // general purpose CPU
        const SAMPLE m_v = m_state[m-1];
        for (int i=m-1; i>0; i--) {
            m_state[i] = (m_state[i-1] + m_v*m_K[i]) % p;
        }
        m_state[0] = (input + m_v*m_K[0]) % p;
        // m_v = m_state[m-1];
#endif
    }
    void back(SAMPLE input=0) {
        const SAMPLE m_v = (m_inv_K0*(m_state[0] - input + p)) % p;
        for (int i=0; i<m-1; i++) {
            m_state[i] = (m_state[i+1] - m_v*m_K[i+1] + p*p) % p;
        }
        m_state[m-1] = m_v;
    }
    void saturate(int q=m) {
        assert(! is_zero<m>(m_state));
        assert(q >= m);
        for (int i=0; i<q; ++i) {
            next();
        }
    }
    bool is_state(STATE st) const {
#ifdef USE_SSE
        bool res = true;
        for (int i=0; i<m; ++i) {
            res &= (m_state[i] == st[i]);
        }
        return res;
#else
        return (st == m_state);
#endif
    }
    auto get_state() const {
        return m_state;
    }
    auto get_cell(int idx) const {
        return m_state[idx];
    }
private:
    alignas(16) STATE m_state {};
    alignas(16) STATE m_K {};
    SAMPLE m_inv_K0 {};
    void m_calc_inv_K0() {
        const auto x = m_K[0];
        assert(x != 0);
        m_inv_K0 = 1;
        while (true) {
            if (((x*m_inv_K0) % p) == 1) {
                break;
            }
            m_inv_K0++;
        }
    }
};

template <int p>
class LFSR_paired_2x4 {
    static_assert(p < 256);
    static_assert(p > 1);
public:
    constexpr LFSR_paired_2x4(u16x8 K): m_K(K) {};
    void set_state(u16x8 state) {
        m_state = state;
    }
    void set_unit_state() {
        m_state = {1, 0, 0, 0, 1, 0, 0, 0};
    }
    void set_K(u16x8 K) {
        m_K = K;
    }
    void next(u16 input=0) {
#ifdef USE_SSE
        __m128i a = _mm_set1_epi16(m_state[3]);
        __m128i b = _mm_set1_epi16(m_state[3] ^ m_state[7]);
        b = _mm_slli_si128(b, 8);
        a = _mm_xor_si128(a, b);
        b = _mm_load_si128((const __m128i*)&m_K[0]);
        __m128i c = _mm_mullo_epi16(a, b);

        const __m128i mask = _mm_set_epi16(-1, -1, -1, 0, -1, -1, -1, 0);
        __m128i inp = _mm_andnot_si128( mask, _mm_set1_epi16(input) );

        __m128i d = _mm_load_si128((const __m128i*)&m_state[0]);
        d = _mm_and_si128(mask, _mm_slli_si128(d, 2));
        d = _mm_add_epi16(c, d);
        d = _mm_add_epi16(inp, d);
        _mm_store_si128((__m128i*)&m_state[0], d);
        for (int i=0; i<8; ++i) {
            m_state[i] %= p;
        }
#else
        u16 m_v3 = m_state[3];
        u16 m_v7 = m_state[7];
        for (int i=7; i>4; i--) {
            m_state[i] = (m_state[i-1] + m_v7*m_K[i]) % p;
            m_state[i-4] = (m_state[i-1-4] + m_v3*m_K[i-4]) % p;
        }
        m_state[0] = (input + m_v3*m_K[0]) % p;
        m_state[4] = (input + m_v7*m_K[4]) % p;
#endif
    }
    void next(u16 inp1, u16 inp2) {
#ifdef USE_SSE
        __m128i a = _mm_set1_epi16(m_state[3]);
        __m128i b = _mm_set1_epi16(m_state[3] ^ m_state[7]);
        b = _mm_slli_si128(b, 8);
        a = _mm_xor_si128(a, b);
        b = _mm_load_si128((const __m128i*)&m_K[0]);
        __m128i c = _mm_mullo_epi16(a, b);

        __m128i mask1 = _mm_slli_si128(_mm_set1_epi16(-1), 2);
        const __m128i mask2 = _mm_set_epi16(-1, -1, -1, 0, -1, -1, -1, -1);
        __m128i inp = _mm_andnot_si128( mask1, _mm_set1_epi16(inp1) );
        inp = _mm_or_si128( inp, _mm_andnot_si128( mask2, _mm_set1_epi16(inp2) ) );

        const __m128i mask = _mm_and_si128(mask1, mask2); // _mm_set_epi16(-1, -1, -1, 0, -1, -1, -1, 0);
        __m128i d = _mm_load_si128((const __m128i*)&m_state[0]);
        d = _mm_and_si128(mask, _mm_slli_si128(d, 2));
        d = _mm_add_epi16(c, d);
        d = _mm_add_epi16(inp, d);
        _mm_store_si128((__m128i*)&m_state[0], d);
        for (int i=0; i<8; ++i) {
            m_state[i] %= p;
        }
#else
        u16 m_v3 = m_state[3];
        u16 m_v7 = m_state[7];
        for (int i=7; i>4; i--) {
            m_state[i] = (m_state[i-1] + m_v7*m_K[i]) % p;
            m_state[i-4] = (m_state[i-1-4] + m_v3*m_K[i-4]) % p;
        }
        m_state[0] = (inp1 + m_v3*m_K[0]) % p;
        m_state[4] = (inp2 + m_v7*m_K[4]) % p;
#endif
    }
    auto get_state() const {
        return m_state;
    }
    bool is_state_low(u16x8 st) const {
        bool res = true;
        for (int i=0; i<4; ++i) {
            res &= (m_state[i] == st[i]);
        }
        return res;
    }
    bool is_state_high(u16x8 st) const {
        bool res = true;
        for (int i=0; i<4; ++i) {
            res &= (m_state[i+4] == st[i+4]);
        }
        return res;
    }
private:
    u16x8 m_state {};
    u16x8 m_K {};
};

}
