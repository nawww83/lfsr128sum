#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#include <clocale>
#endif

#include "lfsr_hash.h"
#include "tests.h"

using namespace std;

constexpr size_t chunkSize = 8*1024*1024;
static vector<uint8_t> buffer(chunkSize, 0);

static lfsr_hash::gens generator;

int main(int argc, char* argv[])
{
#ifdef _WIN32
    std::setlocale(LC_CTYPE, ".UTF8");
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
    if (argc < 2) {
        std::cout << "Передайте путь к файлу при запуске утилиты хэширования:\n";
        std::cout << "> lfsr128sum path_to_file\n";
        std::cout << std::flush;
        //
        lfsr_hash_benchmark();
        // test_lfsr_hash_coverage_1();
        // test_lfsr_hash_coverage_2();
        //
        return 0;
    }
    //
    const std::string& file_name = argv[1];
    std::ifstream fin(file_name, std::ifstream::binary);
    if (!fin) {
        std::cerr << "Ошибка открытия файла: " << file_name << std::endl;
        return 1;
    }
    fin.seekg(0, std::ios::end);
    const auto file_size = fin.tellg();
    assert(file_size > 0);
    fin.seekg(0, std::ios::beg);
    using namespace lfsr_hash;
    constexpr size_t blockSize = 64*1024; // doesn't change it!
    const salt file_salt = {static_cast<int>(file_size % blockSize), static_cast<u16>(file_size), static_cast<u16>(file_size)};
    static_assert((chunkSize % blockSize) == 0);
    u128 hash = {0, 0};
    while (!fin.eof()) {
        fin.read(reinterpret_cast<char*>( buffer.data() ), buffer.size());
        const auto bytesRead = fin.gcount();
        if (bytesRead == chunkSize) {
            constexpr size_t n = chunkSize / blockSize;
            for (size_t i=0; i<n; ++i) {
                u128 inner_hash = hash128<blockSize>(generator, buffer.data() + i*blockSize);
                hash.first ^= inner_hash.first;
                hash.second ^= inner_hash.second;
            }
        } else {
            std::fill(buffer.begin() + (ptrdiff_t)bytesRead,
                      buffer.begin() + (ptrdiff_t)bytesRead + (ptrdiff_t)blockSize, 0); // Make zero padding.
            const size_t n = bytesRead / blockSize;
            const size_t r = bytesRead % blockSize;
            for (size_t i=0; i<n; ++i) {
                u128 inner_hash = hash128<blockSize>(generator, buffer.data() + i*blockSize, file_salt);
                hash.first ^= inner_hash.first;
                hash.second ^= inner_hash.second;
            }
            if (r > 0) {
                u128 inner_hash = hash128<blockSize>(generator, buffer.data() + n*blockSize, file_salt);
                hash.first ^= inner_hash.first;
                hash.second ^= inner_hash.second;
            }
        }
    }
    fin.close();
    std::cout << std::hex << std::setw(16) << std::setfill('0') << hash.first <<
                std::setw(16) << std::setfill('0') << hash.second << std::dec <<
                '\t' << file_name << std::endl;
    return 0;
}
