#pragma once

#include <chrono>

namespace timer_n
{

    class Timer
    {
        std::chrono::time_point<std::chrono::steady_clock> mBegin;

    public:
        Timer() : mBegin(std::chrono::steady_clock::now()) {};
        double elapsed_ns()
        {
            auto mEnd = std::chrono::steady_clock::now();
            return std::chrono::duration_cast<std::chrono::nanoseconds>(mEnd - mBegin).count();
        }
        void reset()
        {
            mBegin = std::chrono::steady_clock::now();
        }
    };

}
