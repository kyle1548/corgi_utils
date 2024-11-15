#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>
#include "../include/fitted_coefficient.hpp"
#include <chrono>


std::pair<double, double> rotatePoint(double x1, double y1, double theta) {
    double x_new = x1 * std::cos(theta) - y1 * std::sin(theta);
    double y_new = x1 * std::sin(theta) + y1 * std::cos(theta);
    return {x_new, y_new};
}

int main() {
    using namespace std::complex_literals; // For 1i
    double l1 = 10;
    std::complex<double> A_l, _1j=1i;
    std::complex<double> theta=90 * M_PI / 180.0;
    
    // double theta = 90 * M_PI / 180.0;

    std::pair<double, double> new_xy;
    double x_new;
    double y_new;
    // Forward kinematics calculations
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i=0; i<10000000; i++){
        A_l = l1 * std::exp(1i * theta);
        // x_new = l1 * std::cos(theta) - 0 * std::sin(theta);
        // y_new = l1 * std::sin(theta) + 0 * std::cos(theta);
    }
    auto end = std::chrono::high_resolution_clock::now();
    // 計算執行時間，並轉換成毫秒
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "time: " << duration.count() << " ms" << std::endl;
    
    std::cout << l1 << A_l << std::endl;
    std::cout << new_xy.first << ", " << new_xy.second << std::endl;

    return 0;
}

