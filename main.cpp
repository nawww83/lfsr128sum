#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#include <clocale>
#endif

#include "lfsr_hash.h"

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
        return -1;
    }
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
    generator.reset();
    const salt file_salt1 = {static_cast<int>(file_size % 31), static_cast<u16>(file_size), static_cast<u16>(3*file_size)};
    generator.add_salt(file_salt1); // To be safe when zero padding made, see below.
    generator.add_salt(S1);
    generator.add_salt(S0);
    generator.add_salt( (file_size % 2) == 0 ? S0 : S4 );
    constexpr size_t blockSize = 64*4; // doesn't change it!
    static_assert((chunkSize % blockSize) == 0);
    while (!fin.eof()) {
        fin.read(reinterpret_cast<char*>( buffer.data() ), buffer.size());
        const auto bytesRead = fin.gcount();
        if (bytesRead == chunkSize) {
            constexpr size_t n = chunkSize / blockSize;
            for (size_t i=0; i<n; ++i) {
                generator.process_input<blockSize>(buffer.data() + i*blockSize);
            }
        } else {
            std::fill(buffer.begin() + (ptrdiff_t)bytesRead,
                      buffer.begin() + (ptrdiff_t)bytesRead + (ptrdiff_t)blockSize, 0); // Make zero padding.
            const size_t n = bytesRead / blockSize;
            const size_t r = bytesRead % blockSize;
            for (size_t i=0; i<n; ++i) {
                generator.process_input<blockSize>(buffer.data() + i*blockSize);
            }
            if (r > 0) {
                generator.process_input<blockSize>(buffer.data() + n*blockSize); // Uses zero padding assumption.
            }
        }
    }
    fin.close();
    generator.add_salt(S0);
    generator.add_salt(S1);
    const salt file_salt2 = {static_cast<int>(file_size % 31), static_cast<u16>(3*file_size), static_cast<u16>(file_size)};
    generator.add_salt(file_salt2); // To be safe when zero padding made, see above.
    generator.add_salt( (file_size % 2) == 0 ? S1 : S0 );
    const u64 h1 = generator.form_hash32();
    generator.add_salt( (file_size % 2) == 0 ? S3 : S2 );
    generator.add_salt( (file_size % 2) == 0 ? S4 : S2 );
    const u64 h2 = generator.form_hash32();
    generator.add_salt( (file_size % 2) == 0 ? S2 : S3 );
    generator.add_salt( (file_size % 2) == 0 ? S3 : S1 );
    const u64 h3 = generator.form_hash32();
    generator.add_salt( (file_size % 2) == 0 ? S4 : S2 );
    generator.add_salt( (file_size % 2) == 0 ? S2 : S0 );
    u64 h4 = generator.form_hash32();
    const u128& hash = {
        h1 | (h2 << 32),
        h3 | (h4 << 32)
    };
    std::cout << std::hex << std::setw(16) << std::setfill('0') << hash.first <<
                std::setw(16) << std::setfill('0') << hash.second << std::dec <<
                '\t' << file_name << std::endl;
    return 0;
}
