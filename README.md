# lfsr128sum

Высокопроизводительный 128-битный хэш-алгоритм на базе LFSR-генераторов.

Проект включает:
- консольную утилиту для вычисления хэша файла;
- библиотеку для использования в C++;
- встроенные тесты корректности и бенчмарк производительности.

## Что это

`lfsr128sum` вычисляет 128-битный хэш файла на основе нескольких LFSR-генераторов с перемешиванием состояний. Основная цель — быстро и стабильно хешировать потоковые данные в памяти и на диске.

## Требования

- C++20
- CMake 3.15+
- процессор с поддержкой AVX2
- Windows/MSVC или Linux/GCC/Clang

## Сборка

```bash
git clone https://github.com/nawww83/lfsr128sum.git
cd lfsr128sum

cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

## Использование

### Вычислить хэш файла

```bash
./lfsr128sum path/to/file
```

Пример вывода:

```text
0123456789abcdef0123456789abcdef  file.bin
```

### Запуск тестов

```bash
./lfsr128sum --test
```

### Запуск бенчмарка

```bash
./lfsr128sum --bench
```

### Версия программы

```bash
./lfsr128sum --version
```

## Использование как библиотеки C++

```cpp
#include "lfsr_hash.h"

int main() {
    lfsr_hash::gens generator;
    generator.reset();

    std::array<std::byte, 1024> data{};
    auto hash = lfsr_hash::hash128(generator, std::span(data));

    // hash.first, hash.second — 64-битные части 128-битного значения
    return 0;
}
```

## Основные возможности

- высокая скорость обработки данных;
- поддержка потокового хэширования;
- встроенные тесты качества и покрытия;
- поддержка Windows и Linux;
- CMake-сборка без внешних зависимостей.

## Лицензия

MIT
