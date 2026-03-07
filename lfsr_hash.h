#pragma once

#include "lfsr.h"
#include "io_utils.h"
#include <span>
#include <utility>


namespace lfsr_hash {

using LFSR251x4 = lfsr8::LFSR_paired_2x4<251>;
using LFSR241x4 = lfsr8::LFSR_paired_2x4<241>;
using STATE = lfsr8::u16x8;
using u16 = lfsr8::u16;
using u32 = lfsr8::u32;
using u64 = lfsr8::u64;
using u128 = std::pair<lfsr8::u64, lfsr8::u64>;

static constexpr STATE K1 = {3, 2, 2, 0, 9, 1, 0, 4};  // p = 251
static constexpr STATE K2 = {7, 0, 4, 3, 7, 3, 2, 0};  // p = 241

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
        g_251x4.power_by(S.q);
        g_241x4.power_by(S.q);
        g_251x4.next(S.s0);
        g_241x4.next(S.s1);
    }

    void process_input(std::span<const std::byte> input);

    auto form_hash32() {
        auto st1 = g_251x4.get_state();
        auto st2 = g_241x4.get_state();
        lfsr8::u32 hash;
        hash  = (st1[0] ^ 1 ^ st1[4]) ^ (st2[0] ^ 3 ^ st2[4]);
        hash <<= 8;
        hash |= (st1[1] ^ 3 ^ st1[5]) ^ (st2[1] ^ 1 ^ st2[5]);
        hash <<= 8;
        hash |= (st1[2] ^ 5 ^ st1[6]) ^ (st2[2] ^ 5 ^ st2[6]);
        hash <<= 8;
        hash |= (st1[3] ^ 3 ^ st1[7]) ^ (st2[3] ^ 3 ^ st2[7]);
        
        return hash;
    }
};

void lfsr_hash::gens::process_input(std::span<const std::byte> input)
{
    const size_t n = input.size();
    if (n < 1) return;

    // 1. Встречное движение (минимум 2 байта)
    for (size_t i = 0; i + 2 <= n; i += 2) 
    {
        u16 left = 0, right = 0;
        io.read_mem(left, input.subspan(i));
        io.read_mem(right, input.subspan(n - 2 - i));
        g_251x4.next(left);
        g_241x4.next(right);
    }

    // 2. Безопасное дополнительное перемешивание
    if (n >= 5) // n-3 >= 2, чтение u16 (2 байта) безопасно
    {
        u16 t1 = 0, t2 = 0;
        io.read_mem(t1, input.subspan(1)); 
        io.read_mem(t2, input.subspan(n - 3)); 

        g_251x4.next(t1 ^ 1);
        g_241x4.next(t2 ^ 1);
        
        g_251x4.back(t1, t2);
        g_241x4.back(~t2, ~t1);
    } 
    else if (n == 4) 
    {
        u16 t1 = 0;
        io.read_mem(t1, input.subspan(1)); 
        g_251x4.next(t1); g_241x4.next(t1);
    }

    // 3. Финальный вклад первого байта
    uint8_t first = static_cast<uint8_t>(input[0]);
    u16 x = static_cast<u16>(first | (first << 8));
    g_251x4.next(x);
    g_241x4.next(x);
}

inline u32 hash32(gens& g, std::span<const std::byte> input)
{
    const auto n = input.size();
    g.add_salt(S1);
    g.add_salt(S0);
    g.process_input(input);
    g.add_salt(S2);
    g.add_salt(S3);
    u32 h = g.form_hash32();
    return h;
}

inline u64 hash64(gens& g, std::span<const std::byte> input)
{
    const auto n = input.size();
    g.add_salt(S1);
    g.process_input(input);
    g.add_salt(S3);
    u64 h1 = g.form_hash32();
    g.add_salt(((n % 2) == 0) ? S1 : S2);
    u64 h2 = g.form_hash32();
    return (h1 << 32) | h2;
}

inline u128 hash128(gens& g, std::span<const std::byte> input)
{
    const auto n = input.size();
    g.add_salt(S1);
    g.add_salt(S0);
    g.add_salt(S4);
    g.process_input(input);
    g.add_salt(S0);
    g.add_salt(S1);
    g.add_salt(S3);

    u64 h1 = g.form_hash32();
    g.add_salt(S3);
    g.add_salt(S4);
    u64 h2 = g.form_hash32();
    g.add_salt(S2);
    g.add_salt(S1);
    u64 h3 = g.form_hash32();
    g.add_salt(S4);
    g.add_salt(S0);
    u64 h4 = g.form_hash32();

    return {
        h1 | (h2 << 32),
        h3 | (h4 << 32)};
}

}