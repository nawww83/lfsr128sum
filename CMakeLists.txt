cmake_minimum_required(VERSION 3.15)
project(lfsr128sum CXX)

# 1. Базовые настройки языка
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 2. ОПРЕДЕЛЕНИЕ ЦЕЛИ (Сначала создаем исполняемый файл)
# Включаем все файлы исходников и заголовков
add_executable(lfsr128sum 
    main.cpp 
    lfsr_hash.h 
    tests.h
)

# 3. НАСТРОЙКА ФЛАГОВ SIMD
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

# 4. ДОПОЛНИТЕЛЬНЫЕ ОПЦИИ СБОРКИ
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Автозапуск тестов после компиляции в режиме отладки
    add_custom_command(TARGET lfsr128sum POST_BUILD
        COMMAND lfsr128sum --test
        COMMENT ">>> Запуск автоматических тестов (Debug) <<<"
    )
else()
    # Максимальная оптимизация для Release
    if(NOT MSVC)
        target_compile_options(lfsr128sum PRIVATE -O3 -ffast-math -ftree-vectorize)
    endif()
endif()

# 5. ТЕСТИРОВАНИЕ (Интеграция с CTest)
enable_testing()
add_test(NAME BaseCoverageTest COMMAND lfsr128sum --test)
add_test(NAME PerformanceBench COMMAND lfsr128sum --bench)

# 6. БИБЛИОТЕКИ
find_package(Threads REQUIRED)
target_link_libraries(lfsr128sum PRIVATE Threads::Threads)
