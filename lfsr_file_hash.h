#pragma once

#include <charconv>
#include <semaphore>
#include <memory>
#include <thread>
#include <span>
#include <filesystem>
#include <string>

#include "lfsr_hash.h"
#include "progress_bar.h"

namespace fs = std::filesystem;

constexpr size_t chunkSize = 8 * 1024 * 1024;
constexpr size_t blockSize = 64 * 1024;


#ifdef _WIN32
#include <malloc.h>
auto aligned_deleter = [](uint8_t *p)
{ _aligned_free(p); };
#define ALLOC_ALIGNED(s) static_cast<uint8_t *>(_aligned_malloc(s, 64))
#else
auto aligned_deleter = [](uint8_t *p)
{ std::free(p); };
#define ALLOC_ALIGNED(s) static_cast<uint8_t *>(std::aligned_alloc(64, s))
#endif

namespace lfsr_file_hash {

using namespace lfsr_hash;

constexpr std::string_view HASH_VERSION_PREFIX = "$lfsr128$v1$";

// Структура для результата парсинга строки хэша
struct ParsedHash {
    std::string_view version;
    std::string_view hex_hash;
    bool is_valid = false;
};

// Извлечение версии алгоритма из строки
inline ParsedHash parse_versioned_hash(std::string_view full_hash_str) {
    ParsedHash result;
    if (full_hash_str.starts_with("$lfsr128$")) {
        size_t version_end = full_hash_str.find('$', 9); // Ищем закрывающий '$' после $lfsr128$
        if (version_end != std::string_view::npos) {
            result.version = full_hash_str.substr(0, version_end + 1);
            result.hex_hash = full_hash_str.substr(version_end + 1);
            result.is_valid = true;
            return result;
        }
    }
    // Обратная совместимость на случай, если префикса нет (чистый hex)
    result.version = "legacy";
    result.hex_hash = full_hash_str;
    result.is_valid = (full_hash_str.size() == 32); // Длина u128 в hex
    return result;
}

// Конвертация строки hex в структуру u128
inline std::optional<u128> hex_to_u128(std::string_view hex_str) {
    if (hex_str.size() != 32) return std::nullopt;

    u128 result;
    auto part1 = hex_str.substr(0, 16);
    auto part2 = hex_str.substr(16, 16);

    // Используем быстрый и безопасный std::from_chars из C++17
    auto res1 = std::from_chars(part1.data(), part1.data() + part1.size(), result.first, 16);
    auto res2 = std::from_chars(part2.data(), part2.data() + part2.size(), result.second, 16);

    if (res1.ec != std::errc{} || res2.ec != std::errc{}) {
        return std::nullopt;
    }
    return result;
}

// Функция принимает путь к файлу, ProgressBar по ссылке и смещение по байтам, и возвращает u128 хэш
inline u128 calculate_file_hash128(const fs::path& p, ProgressBar& bar, uint64_t overall_offset) {
    if (!fs::exists(p) || !fs::is_regular_file(p)) {
        throw std::runtime_error("Файл не найден или недоступен.");
    }

    const uint64_t total_size = fs::file_size(p);
    const salt file_salt = {
        static_cast<int>(total_size % blockSize),
        static_cast<uint16_t>(total_size & 0xFFFF),
        static_cast<uint16_t>((total_size >> 16) ^ (total_size >> 32))
    };

#ifdef _WIN32
    // Флаг 'S' внутри fopen в Windows заставляет файловую систему оптимизировать кэш под последовательный доступ (Sequential)
    FILE *f = fopen(p.string().c_str(), "rbS");
#else
    FILE *f = fopen(p.string().c_str(), "rb");
#endif

    if (!f) throw std::runtime_error("Не удалось открыть файл.");

    auto bufferA_sptr = std::unique_ptr<uint8_t[], decltype(aligned_deleter)>(ALLOC_ALIGNED(chunkSize), aligned_deleter);
    auto bufferB_sptr = std::unique_ptr<uint8_t[], decltype(aligned_deleter)>(ALLOC_ALIGNED(chunkSize), aligned_deleter);
    std::span<uint8_t> bufferA(bufferA_sptr.get(), chunkSize);
    std::span<uint8_t> bufferB(bufferB_sptr.get(), chunkSize);

    size_t bytesInA = 0, bytesInB = 0;
    bool isLastA = false, isLastB = false;

    // Разрешаем диску заполнить до 2 буферов наперед, прежде чем он будет принудительно остановлен
    std::counting_semaphore<2> can_read{2};
    std::counting_semaphore<2> can_process{0};

    u128 total_hash = {0, 0};
    bool done = false;

    // Создаем генератор в куче (heap). Память в куче под x64 всегда выровнена по умолчанию,
    // и объект гарантированно изолирован от стека вызовов функций.
    auto local_generator_ptr = std::make_unique<gens>();
    local_generator_ptr->reset();

    std::thread consumer([&]() {
        uint64_t file_hashed_bytes = 0; // Локальный счетчик обработанных байт
        while (true) {
            can_process.acquire();
            if (done && bytesInA == 0 && bytesInB == 0) break;

            bool processingA = (bytesInA > 0);
            auto& currentBuf = processingA ? bufferA : bufferB;
            size_t currentBytes = processingA ? bytesInA : bytesInB;
            bool currentLast = processingA ? isLastA : isLastB;

            if (currentLast) {
                std::fill(currentBuf.begin() + currentBytes, currentBuf.end(), 0);
                local_generator_ptr->add_salt(file_salt);
            }

            const size_t nBlocks = (currentLast ? (currentBytes + blockSize - 1) / blockSize : chunkSize / blockSize);
            for (size_t i = 0; i < nBlocks; ++i) {
                auto data = currentBuf.subspan(i * blockSize, blockSize);
                u128 res = hash128(*local_generator_ptr, std::as_bytes(data));
                total_hash.first ^= res.first;
                total_hash.second ^= res.second;
            }

            // Буфер успешно обработан процессором. Поток диска уже может читать следующий блок,
            // пока мы спокойно и безболезненно занимаемся выводом на экран!
            file_hashed_bytes += currentBytes;
            bar.update(overall_offset + file_hashed_bytes);

            if (processingA) bytesInA = 0; else bytesInB = 0;
            can_read.release();
            if (currentLast) break;
        }
    });

    done = false;
    uint64_t file_processed = 0; // Сколько байт прочитано конкретно из этого файла

    // Главный поток — Чтение (Producer)
    while (!feof(f)) {
        can_read.acquire();
        bool targetA = (bytesInA == 0);
        auto &targetBuf = targetA ? bufferA : bufferB;

        size_t read = fread(targetBuf.data(), 1, chunkSize, f);
        if (read > 0)
        {
            if (targetA) { bytesInA = read; isLastA = (read < chunkSize); }
            else { bytesInB = read; isLastB = (read < chunkSize); }
            file_processed += read;
            can_process.release();
            if (read < chunkSize) break;
        }
        else {
            can_read.release();
            break;
        }
    }

    // Сигнализируем потоку-потребителю, что новых блоков с диска больше не будет
    done = true;
    can_process.release(); // Подталкиваем consumer, если он заснул в ожидании

    if (consumer.joinable()) consumer.join();
    fclose(f);

    return total_hash;
}

}