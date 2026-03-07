#pragma once

#include <cstdint>
#include <cstring>
#include <span>
#include <bit>
#include <concepts>

namespace io_u {

/**
 * @brief Утилиты для безопасного чтения и записи примитивных типов в байтовые буферы.
 * Поддерживает автоматическую обработку порядка байт (Endianness) и любые контейнеры.
 */
struct io_utils {

    /**
     * @brief Записывает целое число в буфер.
     * При необходимости меняет порядок байт на Little Endian (стандарт для большинства протоколов).
     * 
     * @tparam T Тип записываемого числа (uint16_t, uint32_t и т.д.).
     * @tparam B Тип байта в буфере (uint8_t, std::byte, char).
     * @tparam Extent Размер спана (автоматически выводится для массивов и векторов).
     * 
     * @param x Значение для записи.
     * @param buffer Спан, в который производится запись.
     */
    template<std::integral T, typename B, std::size_t Extent>
        requires (sizeof(B) == 1)
    static void copy_to_mem(T x, std::span<B, Extent> buffer) {
        if (buffer.size() < sizeof(T)) return;

        // В C++20 проверка порядка байт выполняется на этапе компиляции
        if constexpr (std::endian::native == std::endian::big) {
            x = swap_bytes(x);
        }
        std::memcpy(buffer.data(), &x, sizeof(T));
    }

    /**
     * @brief Читает целое число из буфера.
     * Автоматически корректирует порядок байт, если нативная архитектура — Big Endian.
     * 
     * @param x Ссылка на переменную, в которую будет записан результат.
     * @param buffer Спан с исходными данными (может быть константным).
     */
    template<std::integral T, typename B, std::size_t Extent>
        requires (sizeof(B) == 1)
    static void read_mem(T& x, std::span<const B, Extent> buffer) {
        if (buffer.size() < sizeof(T)) return;

        std::memcpy(&x, buffer.data(), sizeof(T));
        if constexpr (std::endian::native == std::endian::big) {
            x = swap_bytes(x);
        }
    }

    /**
     * @brief Вспомогательная перегрузка для удобного вызова с контейнерами (массивы, векторы).
     */
    template<std::integral T, typename Container>
    static void copy_to_mem(T x, Container& c) {
        copy_to_mem(x, std::span{c});
    }

    template<std::integral T, typename Container>
    static void read_mem(T& x, const Container& c) {
        read_mem(x, std::span{c});
    }

private:
    /**
     * @brief Реверс байтов для целых чисел. Аналог std::byteswap из C++23.
     */
    template<std::integral T>
    static constexpr T swap_bytes(T val) noexcept {
        if constexpr (sizeof(T) == 1) return val;
        else if constexpr (sizeof(T) == 2) return __builtin_bswap16(static_cast<uint16_t>(val));
        else if constexpr (sizeof(T) == 4) return __builtin_bswap32(static_cast<uint32_t>(val));
        else if constexpr (sizeof(T) == 8) return __builtin_bswap64(static_cast<uint64_t>(val));
        else {
            T res{};
            for (std::size_t i = 0; i < sizeof(T); ++i) {
                res = (res << 8) | (val & 0xFF);
                val >>= 8;
            }
            return res;
        }
    }
};

} // namespace io_u
