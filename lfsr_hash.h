#pragma once

#include "io_utils.h"
#include "lfsr.h"
#include <utility>

namespace lfsr_hash {

using LFSR251x4 = lfsr8::LFSR_paired_2x4<251>;
using LFSR241x4 = lfsr8::LFSR_paired_2x4<241>;
using STATE = lfsr8::u16x8;
using u16 = lfsr8::u16;
using u32 = lfsr8::u32;
using u64 = lfsr8::u64;
using u128 = std::pair<lfsr8::u64, lfsr8::u64>;

static constexpr STATE K1 = {7, 1, 6, 0, 4, 1, 3, 2};    // p=251
static constexpr STATE K2 = {13, 2, 5, 10, 7, 0, 10, 1}; // p=241

struct salt {
    int q;
    u16 s0;
    u16 s1;
};

static constexpr salt S0 {7, 2, 3};
static constexpr salt S1 {6, 4, 7};
static constexpr salt S2 {31, 8, 11};
static constexpr salt S3 {29, 9, 5};
static constexpr salt S4 {37, 2, 13};

static_assert(S0.q >= 4); // enough saturation, m = 4
static_assert(S1.q >= 4);

static_assert(S2.q >= 6*4); // long distance, m = 4
static_assert(S3.q >= 6*4);
static_assert(S4.q >= 6*4);

io_u::io_utils io;

struct gens {
    LFSR251x4 g_251x4;
    LFSR241x4 g_241x4;
public:
    constexpr gens(): g_251x4(K1),
        g_241x4(K2) {}
    void reset() {
        g_251x4.set_unit_state();
        g_241x4.set_unit_state();
    }
    void add_salt(salt S) {
        for (int i=0; i<S.q; ++i) {
            g_251x4.next(S.s0);
            g_241x4.next(S.s1);
        }
    }

    template <size_t N>
    void process_input(const uint8_t* input) {
        static_assert(N != 0);
        if constexpr (N > 1) {
            for (size_t i=0; i<N/2; ++i) {
                u16 tmp1;
                u16 tmp2;
                io.read_mem_16(tmp1, input + 2*i, sizeof(u16));
                io.read_mem_16(tmp2, input + N - 2 - 2*i, sizeof(u16));
                g_251x4.next(tmp1);
                g_241x4.next(tmp2);
            }
        }
        if constexpr (N > 2) {
            u16 tmp1;
            u16 tmp2;
            io.read_mem_16(tmp1, input + 1, sizeof(u16));
            io.read_mem_16(tmp2, input + N - 3, sizeof(u16));
            g_251x4.next(tmp1);
            g_241x4.next(tmp2);
            g_251x4.next(tmp1);
            g_241x4.next(tmp2);
            g_251x4.next(tmp1);
            g_241x4.next(tmp2);
        }
        u16 x = (u16)input[0] | ((u16)(input[0]) << 8);
        g_251x4.next(x);
        g_241x4.next(x);
    }

    auto form_hash32() {
        auto st1 = g_251x4.get_state();
        auto st2 = g_241x4.get_state();
        lfsr8::u32 hash;
        hash  = ((st1[0] ^ st1[4])) ^ ((st2[0] ^ st2[4]));
        hash <<= 8;
        hash |= ((st1[1] ^ st1[5])) ^ ((st2[1] ^ st2[5]));
        hash <<= 8;
        hash |= ((st1[2] ^ st1[6])) ^ ((st2[2] ^ st2[6]));
        hash <<= 8;
        hash |= ((st1[3] ^ st1[7])) ^ ((st2[3] ^ st2[7]));
        return hash;
    }
};

}
