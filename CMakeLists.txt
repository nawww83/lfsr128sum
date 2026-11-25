cmake_minimum_required(VERSION 3.15)

# 1. Включаем политику для управления рантаймом MSVC
cmake_policy(SET CMP0091 NEW)

project(
    lfsr128sum 
    VERSION 4.0.0
    LANGUAGES CXX
)

# Генерируем реальный .h файл из шаблона во время сборки
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/version.h"
)

# В блоке CPack просто используем эту переменную
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")

# 2. Базовые настройки языка
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 3. ВКЛЮЧАЕМ СТАТИЧЕСКУЮ ЛИНКОВКУ (/MT и /MTd)
# Это избавит от необходимости таскать за собой VCRUNTIME140.dll
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

# 4. Настройка UTF-8 для консоли (чтобы кириллица работала везде)
if(MSVC)
    add_compile_options(/utf-8)
endif()

if(WIN32)
    # Флаги для LLVM-MinGW, чтобы сделать исполняемый файл автономным
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")

    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# 5. ОПРЕДЕЛЕНИЕ ЦЕЛИ (Сначала создаем исполняемый файл)
# Включаем все файлы исходников и заголовков
add_executable(lfsr128sum 
    main.cpp
    progress_bar.h
    timer.h
    lfsr_test_suite.h
    lfsr.h
    lfsr_file_hash.h
    lfsr_hash.h
    hardcores.h
    io_utils.h
)

# Разрешаем компилятору искать заголовок в папке сборки
target_include_directories(lfsr128sum PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

include(CheckCXXSourceRuns)

# Проверка поддержки AVX2 на текущем процессоре
set(AVX2_CODE "
#include <immintrin.h>
int main() {
    __m256i a = _mm256_setzero_si256();
    return 0;
}")

# 1. Временно добавляем флаг для компилятора Clang/GCC
set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
if(MSVC)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} /arch:AVX2")
else()
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -mavx2")
endif()

check_cxx_source_runs("${AVX2_CODE}" HAS_AVX2)

# 2. Возвращаем флаги обратно
set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")

# AVX2 не обязателен: без него собирается эквивалентный SSE4.2-путь
# (результаты хэширования побитно те же, скорость ниже примерно вдвое).
option(LFSR_ARCH_NATIVE "Оптимизировать под CPU сборки (-march=native, AVX2)" ON)

if(HAS_AVX2 AND LFSR_ARCH_NATIVE)
    message(STATUS "SIMD: AVX2, слитный блоковый процессор")
elseif(HAS_AVX2)
    message(STATUS "SIMD: SSE4.2 (LFSR_ARCH_NATIVE=OFF), результаты идентичны AVX2")
else()
    message(WARNING
        "CPU не поддерживает AVX2. Будет собран эквивалентный SSE4.2-путь: "
        "результаты хэширования побитно совпадают; пропускная способность "
        "ниже примерно вдвое (~4 ГБ/с против ~8 ГБ/с у AVX2).")
    message(STATUS "SIMD: SSE4.2")
endif()

# 6. НАСТРОЙКА ФЛАГОВ SIMD
if(MSVC)
    # /Oi включает генерацию внутренних функций (intrinsics).
    # /arch:AVX2 автоматически включает SSE4.2 и ниже, поэтому он нужен
    # только когда AVX2 реально поддерживается.
    target_compile_options(lfsr128sum PRIVATE /Oi)
    if(HAS_AVX2 AND LFSR_ARCH_NATIVE)
        target_compile_options(lfsr128sum PRIVATE /arch:AVX2)
    endif()
    target_compile_definitions(lfsr128sum PRIVATE SIMD_ENABLED=1)
else()
    # -msse4.2 гарантирует SSE 4.2 (нужен для pshufb/blend/mullo_epi32).
    # -march=native добавляет AVX2 и тюнинг под CPU сборки; с
    # -DLFSR_ARCH_NATIVE=OFF получается переносимый SSE4.2-бинарник.
    target_compile_options(lfsr128sum PRIVATE -msse4.2)
    if(HAS_AVX2 AND LFSR_ARCH_NATIVE)
        target_compile_options(lfsr128sum PRIVATE -march=native)
    endif()
    target_compile_definitions(lfsr128sum PRIVATE SIMD_ENABLED=1)
endif()

# 7. ДОПОЛНИТЕЛЬНЫЕ ОПЦИИ СБОРКИ
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Автозапуск тестов после компиляции в режиме отладки
    add_custom_command(TARGET lfsr128sum POST_BUILD
        COMMAND lfsr128sum --test
        COMMENT ">>> Запуск автоматических тестов (Debug) <<<"
    )
else() # Release
    if(MSVC)
        # /O2 - скорость, /Ot - приоритет скорости над размером, /GL - оптимизация всей программы
        target_compile_options(lfsr128sum PRIVATE /O2 /Ot /GL)
        # /LTCG - генерация кода при линковке (нужна для /GL)
        set_property(TARGET lfsr128sum APPEND_STRING PROPERTY LINK_FLAGS " /LTCG")
    else()
        target_compile_options(lfsr128sum PRIVATE -O3 -ffast-math -ftree-vectorize)
    endif()
endif()

# 8. Чтобы консоль Windows понимала UTF-8 через манифест (опционально, но полезно)
if(MSVC)
    set_target_properties(lfsr128sum PROPERTIES LINK_FLAGS "/MANIFEST:EMBED")
endif()

# 9. ТЕСТИРОВАНИЕ (Интеграция с CTest)
enable_testing()
add_test(NAME BaseCoverageTest COMMAND lfsr128sum --test)
add_test(NAME PerformanceBench COMMAND lfsr128sum --bench)

# 10. БИБЛИОТЕКИ
find_package(Threads REQUIRED)
target_link_libraries(lfsr128sum PRIVATE Threads::Threads)

# 11. УПАКОВКА (CPack)
# Указываем, что исполняемый файл идет в корень архива
install(TARGETS lfsr128sum DESTINATION .)

# Добавляем текстовые файлы в корень архива
# Если файлы лежат в подпапках, укажите путь, например: docs/readme.md
install(FILES 
    "${CMAKE_CURRENT_SOURCE_DIR}/README.md" 
    "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE" # или лицензия.txt
    DESTINATION .
)

# Настройки архива
set(CPACK_PACKAGE_NAME "lfsr128sum")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-win64-avx2")
set(CPACK_PACKAGE_VENDOR "Anatoly Novikov")
set(CPACK_GENERATOR "ZIP")
set(CPACK_STRIP_FILES TRUE) 

include(CPack)
