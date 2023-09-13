#include <iostream>
#include <iomanip>
#include <chrono>

class ProgressBar {
private:
    int total;
    int step;
    int barWidth;
    std::chrono::time_point<std::chrono::steady_clock> start_time;

public:
    ProgressBar(int total, int barWidth = 50) : total(total), step(0), barWidth(barWidth) {
        start_time = std::chrono::steady_clock::now();
    }

    void increment() {
        step++;
        int progress = static_cast<int>((static_cast<float>(step) / total) * barWidth);

        std::cout << "[";
        for (int i = 0; i < barWidth; i++) {
            if (i < progress) std::cout << "=";
            else if (i == progress) std::cout << ">";
            else std::cout << " ";
        }

        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_time).count();
        auto estimated_total_time = step == 0 ? 0 : static_cast<int>(elapsed * (static_cast<float>(total) / step));
        auto time_remaining = estimated_total_time - elapsed;

        std::cout << "] "
                  << std::fixed << std::setprecision(2) << (static_cast<float>(step) / total) * 100.0
                  << "% | "
                  << "elapsed: " << elapsed << "s"
                  << " | estimated total: " << estimated_total_time << "s"
                  << " | remaining: " << time_remaining << "s"
                  << "\r";
        std::cout.flush();
    }

    void finish() {
        std::cout << std::endl;
    }
};

