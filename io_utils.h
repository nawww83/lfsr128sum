#pragma once

#include <cstdint>
#include <cassert>
#include <cstring>

namespace io_u {

struct io_utils {

    bool is_little_endian() {
        struct {
            union {
                uint16_t x : 16;
                char b[2];
            } val = {1};
        } t;
        return (t.val.b[0] == 1) && (t.val.b[1] == 0);
    }

    bool is_big_endian() {
        struct {
            union {
                uint16_t x : 16;
                char b[2];
            } val = {1};
        } t;
        return (t.val.b[0] == 0) && (t.val.b[1] == 1);
    }

    void copy_to_mem_16(uint16_t x, uint8_t* buffer, size_t size) {
        static_assert(sizeof(buffer[0]) == 1, "");
        assert(size >= sizeof(uint16_t));
        if (is_little_endian()) {
            std::memcpy(buffer, &x, sizeof(uint16_t));
            return;
        }
        {
            const uint16_t t = x;
            x = ((t >> 8) & 0xff);
            x |= ((t << 8) & 0xff00);
            std::memcpy(buffer, &x, sizeof(uint16_t));
        }
    }

    void copy_to_mem_32(uint32_t x, uint8_t* buffer, size_t size) {
        static_assert(sizeof(buffer[0]) == 1, "");
        assert(size >= sizeof(uint32_t));
        if (is_little_endian()) {
            std::memcpy(buffer, &x, sizeof(uint32_t));
            return;
        }
        {
            const uint32_t t = x;
            x = ((t >> 24) & 0xff);
            x |= ((t >> 8) & 0xff00);
            x |= ((t << 8) & 0xff0000);
            x |= ((t << 24) & 0xff000000);
            std::memcpy(buffer, &x, sizeof(uint32_t));
        }
    }

    void read_mem_16(uint16_t& x, const uint8_t* buffer, size_t size) {
        static_assert(sizeof(buffer[0]) == 1, "");
        assert(size >= sizeof(uint16_t));
        std::memcpy(&x, buffer, sizeof(uint16_t));
        if (! is_little_endian()) {
            const uint16_t t = x;
            x = ((t >> 8) & 0xff);
            x |= ((t << 8) & 0xff00);
        }
    }

    void read_mem_32(uint32_t& x, const uint8_t* buffer, size_t size) {
        static_assert(sizeof(buffer[0]) == 1, "");
        assert(size >= sizeof(uint32_t));
        std::memcpy(&x, buffer, sizeof(uint32_t));
        if (! is_little_endian()) {
            const uint32_t t = x;
            x = ((t >> 24) & 0xff);
            x |= ((t >> 8) & 0xff00);
            x |= ((t << 8) & 0xff0000);
            x |= ((t << 24) & 0xff000000);
        }
    }
}; // struct io_utils

} // io_u
