#pragma once

#ifdef __linux__
#include <pthread.h>
#include <sched.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

namespace hardcores
{

    void bind_to_core(int core_id)
    {
#ifdef __linux__
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(core_id, &cpuset);
        pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#elif defined(_WIN32)
        SetThreadAffinityMask(GetCurrentThread(), 1LL << core_id);
#endif
    }

    template <typename T>
    void doNotOptimizeAway(const T &value)
    {
#if defined(__GNUC__) || defined(__clang__)
        asm volatile("" : : "g"(value) : "memory");
#else
        // Для MSVC — заставляем компилятор думать, что мы модифицируем память
        const volatile T *dummy = &value;
        (void)dummy;
#endif
    }
} // namespace hardcores