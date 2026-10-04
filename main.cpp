#include "progress_bar.h"
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <string_view>

#ifdef _WIN32
#include <windows.h>
#endif

#include "version.h"

#include "lfsr_file_hash.h"
#include "lfsr_test_suite.hpp"

namespace fs = std::filesystem;
using namespace lfsr_hash;

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

#include <cstdlib>

[[maybe_unused]] void print_simd_info()
{
    std::cout << "--- Информация о сборке ---" << std::endl;
    std::cout << "Набор инструкций: " << SIMD_STATUS << std::endl;
    std::cout << "Флаг из CMake:    " << CMAKE_FLAG << std::endl;
    std::cout << "---------------------------" << std::endl;
}

static bool verify_checksum_file(const fs::path& checksum_file_path) {
    fs::path absolute_checksum_path = fs::absolute(checksum_file_path);

    std::ifstream infile(absolute_checksum_path, std::ios::binary); // Открываем как binary, чтобы точно видеть байты
    if (!infile.is_open()) {
        std::cerr << "Ошибка: Не удалось открыть файл контрольных сумм: "
                  << absolute_checksum_path.string() << "\n";
        return false;
    }

    std::string line;
    bool all_files_ok = true;
    size_t line_counter = 0;

    // Счетчики для итогового отчета
    size_t files_checked = 0;
    size_t files_failed = 0;

    // ANSI цвета для вывода
    const std::string_view RESET = "\033[0m";
    const std::string_view RED = "\033[31m";
    const std::string_view GREEN = "\033[32m";
    const std::string_view CYAN = "\033[36m";

    // Добавить в verify_checksum_file сразу после открытия файла:
    char bom[2];
    if (infile.read(bom, 2)) {
        if ((unsigned char)bom[0] == 0xFF && (unsigned char)bom[1] == 0xFE) {
            std::cerr << RED << "Ошибка: Файл хэшей сохранен в кодировке UTF-16 (LE).\n"
                      << "Пожалуйста, пересохраните файл в кодировке UTF-8 или ASCII." << RESET << "\n";
            return false;
        }
    }
    infile.seekg(0); // Сбрасываем указатель чтения в начало, если BOM не найден

    std::cout << CYAN << "* Анализ файла: " << absolute_checksum_path.filename().string() << RESET << "\n";

    while (std::getline(infile, line)) {
        line_counter++;

        // Срезаем Windows/PowerShell BOM (Byte Order Mark) для UTF-8, если он встретился на первой строчке
        if (line_counter == 1 && line.size() >= 3) {
            if ((unsigned char)line[0] == 0xEF && (unsigned char)line[1] == 0xBB && (unsigned char)line[2] == 0xBF) {
                line = line.substr(3);
            }
        }

        // Очищаем от невидимых символов на концах (включая \r и \n)
        while (!line.empty() && (unsigned char)line.back() <= 32) {
            line.pop_back();
        }
        while (!line.empty() && (unsigned char)line.front() <= 32) {
            line.erase(line.begin());
        }

        if (line.empty()) continue;

        // Ищем первый пробел/таб после хэша
        size_t space_pos = line.find_first_of(" \t");
        if (space_pos == std::string::npos) {
            std::cout << RED << "Строка " << line_counter << ": Неверный формат (нет разделителя)" << RESET << "\n";
            all_files_ok = false;
            continue;
        }

        std::string_view raw_hash_part = std::string_view(line).substr(0, space_pos);
        std::string_view filename_part = std::string_view(line).substr(space_pos);

        // Пропускаем все пробелы перед именем файла
        size_t file_pos = filename_part.find_first_not_of(" \t");
        if (file_pos == std::string::npos) {
            std::cout << RED << "Строка " << line_counter << ": Неверный формат (отсутствует имя файла)" << RESET << "\n";
            all_files_ok = false;
            continue;
        }
        filename_part = filename_part.substr(file_pos);

        // 2. Парсим префикс версии хэша
        auto parsed = lfsr_file_hash::parse_versioned_hash(raw_hash_part);
        if (!parsed.is_valid) {
            std::cout << filename_part << ": " << RED << "ИСКАЖЕННАЯ СТРОКА ХЭША" << RESET << "\n";
            all_files_ok = false;
            files_failed++;
            continue;
        }

        // 3. Переводим эталонный хэш из строки в u128
        auto expected_hash_opt = lfsr_file_hash::hex_to_u128(parsed.hex_hash);
        if (!expected_hash_opt) {
            std::cout << filename_part << ": " << RED << "ОШИБКА ФОРМАТА HEX (требуется 32 символа)" << RESET << "\n";
            all_files_ok = false;
            files_failed++;
            continue;
        }
        u128 expected_hash = *expected_hash_opt;

        // 4. Строим абсолютный путь к целевому файлу
        fs::path target_file = fs::absolute(absolute_checksum_path.parent_path() / filename_part);

        if (!fs::exists(target_file)) {
            std::cout << filename_part << ": " << RED << "НЕ НАЙДЕН" << RESET << "\n";
            all_files_ok = false;
            files_failed++;
            continue;
        }

        // 5. Проверяем версию и вычисляем реальный хэш файла
        u128 actual_hash = {0, 0};
        if (parsed.version == lfsr_file_hash::HASH_VERSION_PREFIX || parsed.version == "legacy") {
            try {
                // Создаем временный прогресс-бар для одного конкретного файла.
                // Передаем размер этого файла. Текст лейбла скрываем, так как статус выведется ниже.
                uint64_t f_size = fs::file_size(target_file);
                ProgressBar dummy_bar(f_size, "Checking " + target_file.filename().string());

                // Передаем dummy_bar и смещение 0
                actual_hash = lfsr_file_hash::calculate_file_hash128(target_file, dummy_bar, 0);

                // Стираем строчку прогресс-бара, чтобы она не мешала финальному выводу статуса
                std::cerr << "\r\033[K" << std::flush;
            }
            catch (const std::exception& e) {
                std::cout << filename_part << ": " << RED << "ОШИБКА ЧТЕНИЯ ФАЙЛА" << RESET << "\n";
                all_files_ok = false;
                files_failed++;
                continue;
            }
        } else {
            std::cout << filename_part << ": " << RED << "НЕПОДДЕРЖИВАЕМАЯ ВЕРСИЯ АЛГОРИТМА (" << parsed.version << ")" << RESET << "\n";
            all_files_ok = false;
            files_failed++;
            continue;
        }

        // 6. Сравниваем результаты
        files_checked++;
        if (actual_hash == expected_hash) {
            std::cout << filename_part << ": " << GREEN << "ЦЕЛОСТЕН" << RESET << "\n";
        } else {
            std::cout << filename_part << ": " << RED << "СБОЙ контрольной суммы" << RESET << "\n";
            all_files_ok = false;
            files_failed++;
        }
    }

    // Печатаем красивый финальный отчет
    std::cout << CYAN << "\n=== Итоги проверки ===" << RESET << "\n";
    std::cout << "Всего успешно проверено файлов: " << GREEN << files_checked << RESET << "\n";
    if (files_failed > 0) {
        std::cout << "Файлов с ошибками целостности или формата: " << RED << files_failed << RESET << "\n";
    }

    if (line_counter == 0) {
        std::cout << RED << "Файл контрольных сумм пуст." << RESET << "\n";
        return false;
    }

    return all_files_ok;
}


int main(int argc, char *argv[])
{
    try
    {
#ifdef _WIN32
        SetConsoleOutputCP(CP_UTF8);
        SetConsoleCP(CP_UTF8);
#endif

        if (argc < 2)
        {
            std::cout << "lfsr128sum " << PROJECT_VERSION << "\n";
            std::cout << "Использование:\n"
                      << "  lfsr128sum <файл1> [файл2 ...] [опции]   Вычислить хэш одного или нескольких файлов\n"
                      << "  lfsr128sum -c | --check <файл>     Проверить файлы по контрольным суммам\n"
                      << "Опции:\n"
                      << "  -o, --output <файл>     Записать результат напрямую в указанный файл\n"
                      << "  --test                  Запустить тесты корректности и покрытия\n"
                      << "  --bench                 Запустить бенчмарк производительности\n"
                      << "  --version               Вывести версию программы\n";
            return 0;
        }

        const std::string arg = argv[1];

        if (arg == "--version" || arg == "-v") {
            std::cout << "lfsr128sum version " << PROJECT_VERSION << std::endl;
            return 0;
        }
        if (arg == "--test") {
            LFSRTestSuite suite;
            suite.run_all();
            return 0;
        }
        if (arg == "--bench") {
            LFSRTestSuite suite;
            suite.run_lfsr_benchmark();
            return 0;
        }

        // РЕЖИМ ПРОВЕРКИ ФАЙЛА СУММ
        if (arg == "--check" || arg == "-c") {
            if (argc < 3) {
                std::cerr << "Ошибка: Не указан файл с контрольными суммами.\n";
                return 1;
            }
            fs::path checksum_path(argv[2]);
            bool success = verify_checksum_file(checksum_path);
            return success ? 0 : 1;
        }

        // ====================================================================
        // РЕЖИМ РАСЧЕТА ХЭШЕЙ ДЛЯ МНОЖЕСТВА ФАЙЛОВ (С ИСПРАВЛЕННЫМ PROGRESSBAR)
        // ====================================================================

        std::vector<fs::path> input_files;
        std::string output_file;

        // Разбираем аргументы: ищем файлы и флаг -o / --output
        for (int i = 1; i < argc; ++i) {
            std::string current_arg = argv[i];
            if (current_arg == "-o" || current_arg == "--output") {
                if (i + 1 < argc) { output_file = argv[i + 1]; i++; }
                else { std::cerr << "Ошибка: Не указан файл для записи.\n"; return 1; }
            } else {
                input_files.push_back(fs::path(current_arg));
            }
        }

        if (input_files.empty()) {
            std::cerr << "Ошибка: Не указаны входные файлы.\n"; return 1; }

        // 1. Предварительный проход: отсекаем невалидные файлы и считаем ОБЩИЙ РАЗМЕР
        std::vector<fs::path> valid_files;
        uint64_t total_bytes_to_process = 0;

        for (const auto& file_path : input_files) {
            if (fs::exists(file_path) && fs::is_regular_file(file_path)) {
                valid_files.push_back(file_path);
                total_bytes_to_process += fs::file_size(file_path);
            } else {
                std::cerr << "Предупреждение: Файл не найден или недоступен: " << file_path.string() << "\n";
            }
        }

        if (valid_files.empty()) {
            std::cerr << "Ошибка: Нет доступных файлов для расчета хэша.\n"; return 1; }

        // 2. Создаем ОДИН общий прогресс-бар на весь объем данных
        ProgressBar bar(total_bytes_to_process, "Hashing");
        uint64_t overall_processed_bytes = 0;

        std::stringstream all_results;
        size_t successful_hashes = 0;

        // 3. Последовательно обрабатываем каждый валидный файл
        for (const auto& file_path : valid_files) {
            // Динамически меняем текст лейбла на имя текущего файла
            bar.set_label("Hashing " + file_path.filename().string());

            try {
                // ВАЖНО: Модифицируем calculate_file_hash128, чтобы она принимала
                // ссылку на наш общий bar и счетчик уже обработанных байт
                u128 total_hash = lfsr_file_hash::calculate_file_hash128(file_path, bar, overall_processed_bytes);

                // Накапливаем хэш-сумму
                all_results << lfsr_file_hash::HASH_VERSION_PREFIX
                            << std::hex << std::setw(16) << std::setfill('0') << total_hash.first
                            << std::setw(16) << std::setfill('0') << total_hash.second
                            << std::dec << "  " << file_path.filename().string() << "\n";

                successful_hashes++;
                // Прибавляем размер завершенного файла к общему счетчику прогресса
                overall_processed_bytes += fs::file_size(file_path);
            }
            catch (const std::exception &e) {
                // Если файл упал с ошибкой, всё равно прибавляем его размер, чтобы прогресс-бар не сломался
                overall_processed_bytes += fs::file_size(file_path);
            }
        }

        bar.finish(); // Закрываем общий прогресс-бар (перенос строки в stderr)

        if (successful_hashes == 0) {
            std::cerr << "Ошибка: Не удалось рассчитать хэш ни для одного файла.\n"; return 1; }

        std::string final_output_str = all_results.str();

        // Записываем результат (в файл или консоль)
        if (!output_file.empty()) {
            std::ofstream outfile(output_file, std::ios::out | std::ios::binary);
            if (!outfile.is_open()) {
                std::cerr << "Ошибка записи в файл: " << output_file << "\n"; return 1;
            }
            outfile << final_output_str;
            std::cout << "* Успешно обработано файлов: " << successful_hashes << ". Результаты сохранены в: " << output_file << "\n";
        } else {
            std::cout << final_output_str;
        }

    }
    catch (const std::exception &e)
    {
        std::cerr << "Критическая ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
