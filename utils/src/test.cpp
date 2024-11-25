#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include "../include/fitted_coefficient.hpp"
#include <complex>

#include <chrono>
#include <algorithm>

int main() {
    using namespace std::complex_literals; // For 1i
    double c;
    std::complex<double> a=6i-3.0, b=3i, a1=5i-4.0;
    std::complex<double> theta=90 * M_PI / 180.0;
    std::complex<double> a2;
    // double theta = 90 * M_PI / 180.0;

    std::pair<double, double> new_xy;
    double x = 10.0, y=10.0, z=30.0;

    std::array<std::array<double, 2>, 6> arc_list = {{
        {5.0, 20.0}, 
        {3.0, 15.0},
        {8.0, 10.0},
        {1.0, 25.0},
        {4.0, 30.0},
        {2.0, 18.0}}
    };
    std::array<std::array<double, 2>, 6>::iterator min_arc;
    int min_index;
    double min_value = arc_list[0][0];
    // Forward kinematics calculations
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i=0; i<100000000; i++){
        // c = std::arg((a - a1) / b);
        // a2 = a - a1;
        // c = std::atan2(a2.imag(), a2.real()) - std::atan2(b.imag(), b.real());

        // c = std::atan2(y, x) - std::atan2(-z, 0);
        // c = std::arg((x + 1i*y) / (-1i*z));

        min_arc = std::min_element(arc_list.begin(), arc_list.end(),
                [](const std::array<double, 2>& a, const std::array<double, 2>& b) {
                    return a[0] < b[0];
                });

        // for( int j=1; j<6; j++){
        //     if (arc_list[j][0] < min_value) {
        //     min_value = arc_list[i][0];
        //     min_index = j;
        // }
        // }
    }
    auto end = std::chrono::high_resolution_clock::now();
    // 計算執行時間，並轉換成毫秒
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "time: " << duration.count() << " ms" << std::endl;

    // std::cout << "min_arc = " << min_value << std::endl;

    return 0;
}

