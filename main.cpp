#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <filesystem>
#include <chrono>
#include <thread>
#include <semaphore>

#ifdef _WIN32
#include <windows.h>
#endif

#include "version.h"

#include "lfsr_hash.h"
#include "progress_bar.h"
#include "tests.h"

namespace fs = std::filesystem;
using namespace lfsr_hash;

constexpr size_t chunkSize = 8 * 1024 * 1024;
constexpr size_t blockSize = 64 * 1024;
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

#include <memory>
#include <cstdlib>

#ifdef _WIN32
    #include <malloc.h>
    auto aligned_deleter = [](uint8_t* p) { _aligned_free(p); };
    #define ALLOC_ALIGNED(s) static_cast<uint8_t*>(_aligned_malloc(s, 64))
#else
    auto aligned_deleter = [](uint8_t* p) { std::free(p); };
    #define ALLOC_ALIGNED(s) static_cast<uint8_t*>(std::aligned_alloc(64, s))
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
        // Устанавливаем кодировку UTF-8 для вывода в консоль
        SetConsoleOutputCP(CP_UTF8);
        SetConsoleCP(CP_UTF8);
#endif

        if (argc < 2)
        {
            std::cout << "lfsr128sum " << PROJECT_VERSION << "\n";
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
            std::cout << "lfsr128sum version " << PROJECT_VERSION << std::endl;
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

        // 1. Открываем файл (C-style для скорости)
        FILE* f = fopen(p.string().c_str(), "rb");

        // 2. Выделяем память под системный буфер (8 МБ)
        // Это "пред-буфер", из которого fread будет забирать данные в bufferA/bufferB
        std::vector<char> system_cache(8 * 1024 * 1024); 

        // 3. Устанавливаем полноблочную буферизацию (_IOFBF)
        if (f)
            setvbuf(f, system_cache.data(), _IOFBF, system_cache.size());
        else
            throw std::runtime_error("Не удалось открыть файл.");

        // Два буфера для двойной буферизации
        auto bufferA_sptr = std::unique_ptr<uint8_t[], decltype(aligned_deleter)>(ALLOC_ALIGNED(chunkSize), aligned_deleter);
        auto bufferB_sptr = std::unique_ptr<uint8_t[], decltype(aligned_deleter)>(ALLOC_ALIGNED(chunkSize), aligned_deleter);
        std::span<uint8_t> bufferA(bufferA_sptr.get(), chunkSize);
        std::span<uint8_t> bufferB(bufferB_sptr.get(), chunkSize);

        // Состояние буферов
        size_t bytesInA = 0, bytesInB = 0;
        bool isLastA = false, isLastB = false;

        // Семафоры C++20
        std::binary_semaphore can_read{1};    // Разрешение на чтение (сразу можно)
        std::binary_semaphore can_process{0}; // Разрешение на обработку (ждём данных)

        u128 total_hash = {0, 0};
        bool done = false;
        uint64_t processed = 0;

        // Инициализация генератора перед потоком
        generator.reset();

        // Поток обработки (Consumer)
        std::thread consumer([&]()
                             {
        while (true) {
            can_process.acquire(); // Ждем, пока Producer наполнит какой-то из буферов.

            // Выход, если данных больше нет
            if (done && bytesInA == 0 && bytesInB == 0) break;

            // Определяем, какой из буферов готов (сначала проверяем A, потом B)
            bool processingA = (bytesInA > 0);
            auto& currentBuf = processingA ? bufferA : bufferB;
            size_t currentBytes = processingA ? bytesInA : bytesInB;
            bool currentLast = processingA ? isLastA : isLastB;

            // Логика финализации (Соль и зануление)
            if (currentLast) {
                std::fill(currentBuf.begin() + currentBytes, currentBuf.end(), 0);
                generator.add_salt(file_salt);
            }

            // Цикл хеширования блоками
            const size_t nBlocks = (currentLast ? (currentBytes + blockSize - 1) / blockSize : chunkSize / blockSize);
            for (size_t i = 0; i < nBlocks; ++i) {
                auto data = currentBuf.subspan(i * blockSize, blockSize);
                u128 res = hash128(generator, std::as_bytes(data));
                total_hash.first ^= res.first;
                total_hash.second ^= res.second;
            }

            // Освобождаем буфер, помечая его пустым
            if (processingA) bytesInA = 0; else bytesInB = 0;

            can_read.release(); // Сигнализируем Producer'у, что место освободилось
            if (currentLast) break;
        } });

        ProgressBar bar(total_size, "Hashing ");
        // Главный поток — Чтение (Producer)
        while (!feof(f))
        {
            can_read.acquire(); // Ждем свободный буфер (A или B)

            // Выбираем пустой буфер
            bool targetA = (bytesInA == 0);
            auto &targetBuf = targetA ? bufferA : bufferB;

            size_t read = fread(targetBuf.data(), 1, chunkSize, f);
            if (read > 0)
            {
                if (targetA)
                {
                    bytesInA = read;
                    isLastA = (read < chunkSize);
                }
                else
                {
                    bytesInB = read;
                    isLastB = (read < chunkSize);
                }

                processed += read;
                bar.update(processed);

                can_process.release(); // Сигнализируем Consumer'у, что данные готовы
                if (read < chunkSize)
                    break; // Конец файла
            }
            else
            {
                can_read.release(); // Если ничего не прочитали, возвращаем семафор
                break;
            }
        }

        done = true;
        if (consumer.joinable())
            consumer.join();

        bar.finish();

        fclose(f);

        // Итоговый вывод
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
