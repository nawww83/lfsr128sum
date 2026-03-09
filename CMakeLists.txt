cmake_minimum_required(VERSION 3.15)

# 1. Включаем политику для управления рантаймом MSVC
cmake_policy(SET CMP0091 NEW)

project(
    lfsr128sum 
    VERSION 2.0.0 
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

# 5. ОПРЕДЕЛЕНИЕ ЦЕЛИ (Сначала создаем исполняемый файл)
# Включаем все файлы исходников и заголовков
add_executable(lfsr128sum 
    main.cpp
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

check_cxx_source_runs("${AVX2_CODE}" HAS_AVX2)

if(NOT HAS_AVX2)
    message(FATAL_ERROR "CPU не поддерживает AVX2. Сборка невозможна, так как проект оптимизирован под этот набор инструкций.")
endif()

# 6. НАСТРОЙКА ФЛАГОВ SIMD
if(MSVC)
    # /arch:AVX2 автоматически включает SSE4.2 и ниже. 
    # /Oi включает генерацию внутренних функций (intrinsics).
    target_compile_options(lfsr128sum PRIVATE /arch:AVX2 /Oi)
    target_compile_definitions(lfsr128sum PRIVATE SIMD_ENABLED=1)
else()
    # -march=native заставляет GCC/Clang использовать все инструкции текущего CPU.
    # -msse4.2 гарантирует поддержку SSE 4.2 даже если native не сработал.
    target_compile_options(lfsr128sum PRIVATE -march=native -msse4.2)
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
