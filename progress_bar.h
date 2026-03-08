#pragma once
#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif
#endif

class ProgressBar {
private:
    size_t total_size;
    std::string label;
    int last_pct;
    std::chrono::steady_clock::time_point start_time;
    double smoothed_speed = 0; // Для плавности
    const double alpha = 0.1;   // Коэффициент сглаживания (0.1 - очень плавно, 0.5 - быстро)

public:
    ProgressBar(size_t total, std::string label )
        : total_size(total), label(label), last_pct(-1) {
        
        #ifdef _WIN32
        // Активируем цвета для Windows
        HANDLE hOut = GetStdHandle(STD_ERROR_HANDLE);
        DWORD dwMode = 0;
        if (hOut != INVALID_HANDLE_VALUE && GetConsoleMode(hOut, &dwMode)) {
            SetConsoleMode(hOut, dwMode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
        }
        #endif
        
        start_time = std::chrono::steady_clock::now();
    }

    void update(size_t processed) {
        if (total_size == 0) return;
        int pct = static_cast<int>((processed * 100) / total_size);

        if (pct != last_pct) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - start_time).count();
            
            // Текущая средняя скорость
            double current_speed = (elapsed > 0) ? (processed / elapsed) : 0;
            
            // Формула EMA: сглаживаем текущее значение на основе предыдущего
            if (smoothed_speed == 0) smoothed_speed = current_speed;
            else smoothed_speed = (alpha * current_speed) + (1.0 - alpha) * smoothed_speed;

            double eta = (smoothed_speed > 0) ? (static_cast<double>(total_size - processed) / smoothed_speed) : 0;

            // Собираем всё в одну строку-буфер
            std::stringstream ss;
            ss << "\r\033[K" << label << ": [";
            
            int barWidth = 25;
            int pos = (barWidth * pct) / 100;

            ss << "\033[32m"; // Green
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos) ss << "█";
                else if (i == pos) ss << "▓";
                else ss << "░";
            }
            ss << "\033[0m " << pct << "% ";

            ss << "\033[36m| " << std::fixed << std::setprecision(1) 
               << (smoothed_speed / 1024 / 1024) << " MB/s "
               << "| ETA: " << static_cast<int>(eta) << "s\033[0m";

            // Однократный вывод всей строки убирает фликкер
            std::string output = ss.str();
            std::cerr.write(output.c_str(), output.size());
            std::cerr.flush();

            last_pct = pct;
        }
    }
    void finish() {
        std::cerr << std::endl;
    }
};
