#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <stdexcept>


#include "../include/fitted_coefficient.hpp"
#include <chrono>
#include <Eigen/Dense>


class LegModel {
    public:
        // Constructor
        LegModel(bool sim = true);

        // Forward kinematics
        void forward(double theta, double beta);

        // Inverse kinematics (partial implementation)
        std::pair<double, double> inverse(const std::pair<double, double>& pos, const std::string& joint = "G");

        // Move function (partial implementation)
        std::pair<double, double> move(double theta, double beta, const std::pair<double, double>& move_vec, bool contact_upper = true);


        // Joint Positions
        std::array<double, 2> A_l, A_r, B_l, B_r, C_l, C_r, D_l, D_r, E, F_l, F_r, G, H_l, H_r, U_l, U_r, L_l, L_r;

    private:
        // Joint Positions in complex
        std::complex<double> A_l_c, A_r_c, B_l_c, B_r_c, C_l_c, C_r_c, D_l_c, D_r_c, E_c, F_l_c, F_r_c, G_c, H_l_c, H_r_c, U_l_c, U_r_c, L_l_c, L_r_c;

        // Constants
        double max_theta;
        double min_theta;
        double theta0;
        double beta0;
        double R; // Wheel radius
        double r; // Tire radius
        double radius;
        double arc_HF;
        double arc_BC;
        double l1, l2, l3, l4, l5, l6, l7, l8;
        double l_AE, l_BF, l_BH, ang_UBC, ang_LFG;

        // Intermediate variables
        double l_BD;
        double ang_OEA;
        double ang_BCF;
        double ang_DBC;
        double ang_OGF;

        // Current theta and beta
        double theta;
        double beta;

        // Helper functions
        void calculate();
        void rotate();
        void symmetry();
        void to_vector();
};

LegModel::LegModel(bool sim) {
    // Constants initialization
    max_theta = M_PI * 160.0 / 180.0;
    min_theta = M_PI * 17.0 / 180.0;
    theta0 = M_PI * 17.0 / 180.0;
    beta0 = M_PI * 90.0 / 180.0;

    // Wheel radius
    R = 0.1; // 10 cm
    if (sim) {
        r = 0.0125; // No tire
    } else {
        r = 0.019;  // With tire
    }
    radius = R + r;

    // Linkage parameters
    arc_HF = M_PI * 130.0 / 180.0;
    arc_BC = M_PI * 101.0 / 180.0;
    l1 = 0.8 * R;
    l2 = R - l1;
    l3 = 2.0 * R * sin(arc_BC / 2.0);
    l4 = 0.882966335 * R;
    l5 = 0.9 * R;
    l6 = 0.4 * R;
    l7 = 2.0 * R * sin((arc_HF - arc_BC - theta0) / 2.0);
    l8 = 2.0 * R * sin((M_PI - arc_HF) / 2.0);

    // Useful parameters
    l_AE = l5 + l6;
    l_BF = 2.0 * R * sin((arc_HF - theta0) / 2.0);
    l_BH = 2.0 * R * sin(theta0 / 2.0);
    ang_UBC = (M_PI - arc_BC) / 2.0;
    ang_LFG = (M_PI - (M_PI - arc_HF)) / 2.0;

    // Initialize positions
    forward(theta0, 0.0);
}

void LegModel::forward(double theta_in, double beta_in) {
    theta = theta_in;
    beta = beta_in;

    // Limit theta
    if (theta > max_theta) {
        theta = max_theta;
        std::cout << "Theta exceeds upper limit. Set to max_theta.\n";
    }
    if (theta < min_theta) {
        theta = min_theta;
        std::cout << "Theta below lower limit. Set to min_theta.\n";
    }

    // Calculate positions
    calculate();
    rotate();
    to_vector();
}

void LegModel::calculate() {
    using namespace std::complex_literals; // For 1i
    
    // Forward kinematics calculations
    A_l_c = l1 * std::exp(1i * theta);
    B_l_c = R * std::exp(1i * theta);
    ang_OEA = std::asin((l1 / l_AE) * sin(theta));
    E_c = l1 * cos(theta) - l_AE * cos(ang_OEA);
    D_l_c = E_c + l6 * std::exp(1i * ang_OEA);
    l_BD = std::abs(D_l_c - B_l_c);
    ang_DBC = std::acos((l_BD * l_BD + l3 * l3 - l4 * l4) / (2.0 * l_BD * l3));
    C_l_c = B_l_c + (D_l_c - B_l_c) * std::exp(-1i * ang_DBC) * (l3 / l_BD);
    ang_BCF = std::acos((l3 * l3 + l7 * l7 - l_BF * l_BF) / (2.0 * l3 * l7));
    F_l_c = C_l_c + (B_l_c - C_l_c) * std::exp(-1i * ang_BCF) * (l7 / l3);
    ang_OGF = std::asin(std::abs(F_l_c.imag()) / l8);
    G_c = F_l_c - l8 * std::exp(1i * ang_OGF);
    U_l_c = B_l_c + (C_l_c - B_l_c) * std::exp(1i * ang_UBC) * (R / l3);
    L_l_c = F_l_c + (G_c - F_l_c) * std::exp(1i * ang_LFG) * (R / l8);
    H_l_c = U_l_c + (B_l_c - U_l_c) * std::exp(-1i * theta0);

    symmetry();
}

void LegModel::symmetry() {
    // Symmetric positions
    A_r_c = std::conj(A_l_c);
    B_r_c = std::conj(B_l_c);
    C_r_c = std::conj(C_l_c);
    D_r_c = std::conj(D_l_c);
    F_r_c = std::conj(F_l_c);
    H_r_c = std::conj(H_l_c);
    U_r_c = std::conj(U_l_c);
    L_r_c = std::conj(L_l_c);
}

void LegModel::rotate() {
    using namespace std::complex_literals;
    std::complex<double> rot_ang = std::exp(1i * (beta + beta0));

    // Rotate positions
    A_l_c *= rot_ang;
    A_r_c *= rot_ang;
    B_l_c *= rot_ang;
    B_r_c *= rot_ang;
    C_l_c *= rot_ang;
    C_r_c *= rot_ang;
    D_l_c *= rot_ang;
    D_r_c *= rot_ang;
    E_c   *= rot_ang;
    F_l_c *= rot_ang;
    F_r_c *= rot_ang;
    G_c   *= rot_ang;
    H_l_c *= rot_ang;
    H_r_c *= rot_ang;
    U_l_c *= rot_ang;
    U_r_c *= rot_ang;
    L_l_c *= rot_ang;
    L_r_c *= rot_ang;
}

void LegModel::to_vector() {
    A_l = {A_l_c.real(), A_l_c.imag()};
    A_r = {A_r_c.real(), A_r_c.imag()};
    B_l = {B_l_c.real(), B_l_c.imag()};
    B_r = {B_r_c.real(), B_r_c.imag()};
    C_l = {C_l_c.real(), C_l_c.imag()};
    C_r = {C_r_c.real(), C_r_c.imag()};
    D_l = {D_l_c.real(), D_l_c.imag()};
    D_r = {D_r_c.real(), D_r_c.imag()};
    E   = {E_c.real()  , E_c.imag()};
    F_l = {F_l_c.real(), F_l_c.imag()};
    F_r = {F_r_c.real(), F_r_c.imag()};
    G   = {G_c.real()  , G_c.imag()};
    H_l = {H_l_c.real(), H_l_c.imag()};
    H_r = {H_r_c.real(), H_r_c.imag()};
    U_l = {U_l_c.real(), U_l_c.imag()};
    U_r = {U_r_c.real(), U_r_c.imag()};
    L_l = {L_l_c.real(), L_l_c.imag()};
    L_r = {L_r_c.real(), L_r_c.imag()};
}

// Note: The inverse and move functions require root-finding and numerical methods that are complex to implement.
// For a complete implementation, you would need to use numerical libraries like Eigen, Ceres Solver, or write custom solvers.

int main() {
    LegModel legmodel(true);

    // Forward kinematics example
    std::cout << "****************************************\n";
    std::cout << "****** Forward kinematics example ******\n";
    std::cout << "****************************************\n";

    // Single input
    std::cout << "==========Single Input==========\n";
    double theta = M_PI * 130.0 / 180.0;
    double beta = M_PI;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<1000000; i++){
        legmodel.forward(theta, beta);
    }
    auto end = std::chrono::high_resolution_clock::now();
    // 計算執行時間，並轉換成毫秒
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "time: " << duration.count() << " ms" << std::endl;


    std::cout << "Output G with single value input: (" << legmodel.G[0] << ", " << legmodel.G[1] << ")\n";

    // Note: The contact_map function and other advanced features are not fully implemented in this example.


    Eigen::VectorXd a(5);
    Eigen::VectorXd b(5);
    Eigen::VectorXd result(5);

    // 初始化向量
    a << 1.0, 2.0, 3.0, 4.0, 5.0;

    b << 5.0, 4.0, 3.0, 2.0, 1.0;
    a << 5.0, 2.0, 3.0, 4.0, 5.0;
    std::cout << a ;
    return 0;
}