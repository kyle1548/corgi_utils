#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>


#include "../include/fitted_coefficient.hpp"
#include <chrono>



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

    private:
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

        // Positions
        std::complex<double> A_l, A_r, B_l, B_r, C_l, C_r, D_l, D_r, E, F_l, F_r, G, H_l, H_r, U_l, U_r, L_l, L_r;

        // Current theta and beta
        double theta;
        double beta;

        // Helper functions
        void calculate();
        void rotate();
        void symmetry();
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
    }

    void LegModel::calculate() {
        using namespace std::complex_literals; // For 1i
        
        auto start = std::chrono::high_resolution_clock::now();

        // Forward kinematics calculations
        A_l = l1 * std::exp(1i * theta);
        B_l = R * std::exp(1i * theta);
        ang_OEA = std::asin((l1 / l_AE) * sin(theta));
        E = l1 * cos(theta) - l_AE * cos(ang_OEA);
        D_l = E + l6 * std::exp(1i * ang_OEA);
        l_BD = std::abs(D_l - B_l);
        ang_DBC = std::acos((l_BD * l_BD + l3 * l3 - l4 * l4) / (2.0 * l_BD * l3));
        C_l = B_l + (D_l - B_l) * std::exp(-1i * ang_DBC) * (l3 / l_BD);
        ang_BCF = std::acos((l3 * l3 + l7 * l7 - l_BF * l_BF) / (2.0 * l3 * l7));
        F_l = C_l + (B_l - C_l) * std::exp(-1i * ang_BCF) * (l7 / l3);
        ang_OGF = std::asin(std::abs(F_l.imag()) / l8);
        G = F_l - l8 * std::exp(1i * ang_OGF);
        U_l = B_l + (C_l - B_l) * std::exp(1i * ang_UBC) * (R / l3);
        L_l = F_l + (G - F_l) * std::exp(1i * ang_LFG) * (R / l8);
        H_l = U_l + (B_l - U_l) * std::exp(-1i * theta0);


        auto end = std::chrono::high_resolution_clock::now();
        // 計算執行時間，並轉換成毫秒
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "執行時間: " << duration.count() << " 毫秒" << std::endl;

        symmetry();
    }

    void LegModel::symmetry() {
        // Symmetric positions
        A_r = std::conj(A_l);
        B_r = std::conj(B_l);
        C_r = std::conj(C_l);
        D_r = std::conj(D_l);
        F_r = std::conj(F_l);
        H_r = std::conj(H_l);
        U_r = std::conj(U_l);
        L_r = std::conj(L_l);
    }

    void LegModel::rotate() {
        using namespace std::complex_literals;
        std::complex<double> rot_ang = std::exp(1i * (beta + beta0));

        // Rotate positions
        A_l *= rot_ang;
        A_r *= rot_ang;
        B_l *= rot_ang;
        B_r *= rot_ang;
        C_l *= rot_ang;
        C_r *= rot_ang;
        D_l *= rot_ang;
        D_r *= rot_ang;
        E *= rot_ang;
        F_l *= rot_ang;
        F_r *= rot_ang;
        G *= rot_ang;
        H_l *= rot_ang;
        H_r *= rot_ang;
        U_l *= rot_ang;
        U_r *= rot_ang;
        L_l *= rot_ang;
        L_r *= rot_ang;
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
    double beta = 0.0;
    legmodel.forward(theta, beta);
    // std::cout << "Output G with single value input: (" << legmodel.G.real() << ", " << legmodel.G.imag() << ")\n";

    // Note: The contact_map function and other advanced features are not fully implemented in this example.

    return 0;
}