#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <filesystem>
#include <chrono>
#ifdef _WIN32
#include <windows.h>
#endif

#include "lfsr_hash.h"
#include "tests.h"

static constexpr auto VERSION = "v2.1-simd";

namespace fs = std::filesystem;
using namespace lfsr_hash;

constexpr size_t chunkSize = 8 * 1024 * 1024;
constexpr size_t blockSize = 64 * 1024;
static std::vector<uint8_t> buffer(chunkSize);
static gens generator;

// Проверяем все возможные макросы компиляторов
#if defined(__AVX2__)
#define SIMD_STATUS "AVX2"
#elif defined(__AVX__)
#define SIMD_STATUS "AVX"
#elif defined(__SSE4_2__)
#define SIMD_STATUS "SSE4.2"
#elif defined(__SSE4_1__)
#define SIMD_STATUS "SSE4.1"
#elif defined(_M_AMD64) || defined(_M_X64) || defined(__x86_64__)
#define SIMD_STATUS "x86_64 (Base SSE2)"
#else
#define SIMD_STATUS "НЕ ОПРЕДЕЛЕН (Generic C++)"
#endif

#ifdef SIMD_ENABLED
#define CMAKE_FLAG "ОК"
#else
#define CMAKE_FLAG "ОТСУТСТВУЕТ"
#endif

[[maybe_unused]] void print_simd_info()
{
    std::cout << "--- Информация о сборке ---" << std::endl;
    std::cout << "Набор инструкций: " << SIMD_STATUS << std::endl;
    std::cout << "Флаг из CMake:    " << CMAKE_FLAG << std::endl;
    std::cout << "---------------------------" << std::endl;
}

int main(int argc, char *argv[])
{
    // print_simd_info();
    try
    {
#ifdef _WIN32
        SetConsoleOutputCP(CP_UTF8);
#endif

        if (argc < 2)
        {
            std::cout << "lfsr128sum " << VERSION << "\n"; // <-- Вывод здесь
            std::cout << "Использование: lfsr128sum <путь_к_файлу> [опции]\n\n"
                      << "Опции:\n"
                      << "  --test    Запустить тесты корректности и покрытия\n"
                      << "  --bench   Запустить бенчмарк производительности\n"
                      << "  --version Вывести версию программы\n";
            return 0;
        }

        const std::string arg = argv[1];

        if (arg == "--version" || arg == "-v")
        {
            std::cout << "lfsr128sum version " << VERSION << std::endl;
            return 0;
        }

        if (arg == "--test")
        {
            check_hash();
            test_simd_consistency();
            test_block_simd_consistency();
            test_lfsr_hash_coverage_1();
            test_lfsr_hash_coverage_2();
            return 0;
        }
        if (arg == "--bench")
        {
            lfsr_hash_benchmark();
            return 0;
        }

        fs::path p(arg);
        if (!fs::exists(p) || !fs::is_regular_file(p))
        {
            std::cerr << "Ошибка: Файл не найден или недоступен.\n";
            return 1;
        }

        const uint64_t total_size = fs::file_size(p);
        const salt file_salt = {
            static_cast<int>(total_size % blockSize),
            static_cast<uint16_t>(total_size & 0xFFFF),
            static_cast<uint16_t>((total_size >> 16) ^ (total_size >> 32))};

        std::ifstream fin(p, std::ifstream::binary);
        if (!fin)
            throw std::runtime_error("Не удалось открыть файл.");

        u128 total_hash = {0, 0};
        uint64_t processed = 0;
        int last_pct = -1;
        generator.reset();

        while (fin.read(reinterpret_cast<char *>(buffer.data()), buffer.size()) || fin.gcount() > 0)
        {
            const size_t bytesRead = static_cast<size_t>(fin.gcount());
            const bool isLast = (bytesRead < chunkSize);

            if (isLast)
            {
                std::fill(buffer.begin() + bytesRead, buffer.end(), 0);
                generator.add_salt(file_salt);
            }

            const size_t nBlocks = (isLast ? (bytesRead + blockSize - 1) / blockSize : chunkSize / blockSize);
            for (size_t i = 0; i < nBlocks; ++i)
            {
                auto data = std::span(buffer).subspan(i * blockSize, blockSize);
                u128 res = hash128(generator, std::as_bytes(data));
                total_hash.first ^= res.first;
                total_hash.second ^= res.second;
            }

            processed += bytesRead;
            int pct = static_cast<int>((processed * 100) / total_size);
            if (pct != last_pct && total_size > chunkSize)
            {
                std::cerr << "\rProgress: " << pct << "%" << std::flush;
                last_pct = pct;
            }
        }
        if (total_size > chunkSize)
            std::cerr << "\r" << std::string(20, ' ') << "\r";

        std::cout << std::hex << std::setw(16) << std::setfill('0') << total_hash.first
                  << std::setw(16) << std::setfill('0') << total_hash.second
                  << std::dec << "  " << p.filename().string() << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Критическая ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
