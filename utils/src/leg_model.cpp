#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <stdexcept>


#include "../include/fitted_coefficient.hpp"
#include <chrono>
// #include <Eigen/Dense>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>

class LegModel {
    public:
        // Constructor
        LegModel(bool sim = true);

        // Forward kinematics
        void forward(double theta, double beta, bool vector = true);

        // Inverse kinematics (partial implementation)
        void inverse(const double pos[2], const std::string &joint = "G", bool forward = true);

        void contact_map(double theta_in, double beta_in, double slope = 0);


        // Move function (partial implementation)
        std::pair<double, double> move(double theta, double beta, const std::pair<double, double>& move_vec, bool contact_upper = true);


        // Joint Positions
        std::array<double, 2> A_l, A_r, B_l, B_r, C_l, C_r, D_l, D_r, E, F_l, F_r, G, H_l, H_r, U_l, U_r, L_l, L_r;

        // Current theta and beta
        double theta;
        double beta;

        // Contact map variable
        int rim = 3;    // 1 -> 2 -> 3 -> 4 -> 5 -> 0: 
                        // U_l -> L_l -> G -> L_r -> U_r -> None
        double alpha;

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

        // Helper functions
        void calculate();
        void rotate();
        void symmetry();
        void to_vector();
        std::array<double, 2> arc_min(const std::complex<double>& p1, const std::complex<double>& p2, const std::complex<double>& O, const std::string& rim);
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
    this->forward(theta0, 0.0);
}

void LegModel::forward(double theta_in, double beta_in, bool vector) {
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
    this->calculate();
    this->rotate();
    if (vector) {
        this->to_vector();
    }
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
    this->symmetry();
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
void LegModel::inverse(const double pos[2], const std::string &joint, bool forward) {
        using namespace std::complex_literals;
        double abs_pos = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
        if (joint == "G"){
            theta = inv_G_dist_poly(abs_pos);
            beta = std::atan2(pos[1], pos[0]) - std::atan2(-abs_pos, 0);    // atan2(y, x)
        } else if (joint == "Ul" || joint == "Ur"){
            theta = inv_U_dist_poly(abs_pos);
            double U_x_beta0, U_y_beta0;
            if (joint == "Ul"){
                U_x_beta0 = U_l_poly[0](theta);
                U_y_beta0 = U_l_poly[1](theta);
            } else {    // Ur
                U_x_beta0 = U_r_poly[0](theta);
                U_y_beta0 = U_r_poly[1](theta);    
            }//end if else
            beta = std::atan2(pos[1], pos[0]) - std::atan2(U_y_beta0, U_x_beta0);    // atan2(y, x)
        } else if (joint == "Ll" || joint == "Lr"){
            theta = inv_L_dist_poly(abs_pos);
            double L_x_beta0, L_y_beta0;
            if (joint == "Ll"){
                L_x_beta0 = L_l_poly[0](theta);
                L_y_beta0 = L_l_poly[1](theta);
            } else {    // Lr
                L_x_beta0 = L_r_poly[0](theta);
                L_y_beta0 = L_r_poly[1](theta);  
            }//end if else            
            beta = std::atan2(pos[1], pos[0]) - std::atan2(L_y_beta0, L_x_beta0);    // atan2(y, x)
        } else {
            throw std::runtime_error("joint needs to be 'G', 'Ul', 'Ur', 'Ll', or 'Lr'.");
        }//end if else
        if (forward) {
            this->forward(theta, beta);
        }
}//end inverse


void LegModel::contact_map(double theta_in, double beta_in, double slope) {
        using namespace std::complex_literals;
        double beta_adjusted = beta_in - slope;

        this->forward(theta, beta_adjusted, false);

        std::complex<double> G_l_tmp = (G_c - L_l_c) / R * radius + L_l_c;
        std::complex<double> G_r_tmp = (G_c - L_r_c) / R * radius + L_r_c;
        std::complex<double> H_l_tmp = (H_l_c - U_l_c) / R * radius + U_l_c;
        std::complex<double> H_r_tmp = (H_r_c - U_r_c) / R * radius + U_r_c;
        std::complex<double> F_l_tmp = (F_l_c - U_l_c) / R * radius + U_l_c;
        std::complex<double> F_r_tmp = (F_r_c - U_r_c) / R * radius + U_r_c;

        std::array<std::array<double, 2>, 6> arc_list = {
            this->arc_min(H_l_tmp, F_l_tmp, U_l_c, "left upper"),
            this->arc_min(F_l_tmp, G_l_tmp, L_l_c, "left lower"),
            this->arc_min(G_l_tmp, G_r_tmp, G_c, "G"),
            this->arc_min(G_r_tmp, F_r_tmp, L_r_c, "right lower"),
            this->arc_min(F_r_tmp, H_r_tmp, U_r_c, "right upper"),
            {0.0, 0.0}
        };

        double min_value = arc_list[0][0];
        int min_index = 0;
        for(int i=1; i<6; i++){
            if (arc_list[i][0] < min_value) {
                min_value = arc_list[i][0];
                min_index = i;
            }//end if
        }//end for

        rim = min_index==5? 0 : min_index+1;
        alpha = arc_list[min_index][1];
}//end contact_map

std::array<double, 2> LegModel::arc_min(const std::complex<double>& p1, const std::complex<double>& p2, const std::complex<double>& O, const std::string& rim) {
        using namespace std::complex_literals;
        double lowest_point = 0.0;
        double alpha = 0.0;
        double bias_alpha = 0.0;

        if (rim == "left upper") {
            // bias_alpha = -M_PI;
        } else if (rim == "left lower") {
            // bias_alpha = -M_PI / 3.6; // -50 degrees
        } else if (rim == "G") {
            std::complex<double> direction_G = p1 + p2;
            bias_alpha = std::arg((p1 - O) / direction_G);
        } else if (rim == "right lower") {
            // bias_alpha = 0.0;
        } else if (rim == "right upper") {
            // bias_alpha = M_PI / 3.6; // 50 degrees
        }//end if else

        double cal_err = 1e-9;
        bool in_range = ((p2 - O).real() >= -cal_err) && ((p1 - O).real() <= cal_err);

        if (in_range) {
            lowest_point = O.imag() - radius;
            alpha = std::arg(-1i / (p1 - O));
        } else {
            std::complex<double> smaller = (p1.imag() < p2.imag()) ? p1 : p2;
            lowest_point = 1.0; // Set to a large value if not normal contact
            alpha = std::arg((smaller - O) / (p1 - O));
        }//end if else

        return {lowest_point, alpha + bias_alpha};
}//end arc_min


int main() {
    LegModel legmodel(true);

    // Forward kinematics example
    std::cout << "****************************************\n";
    std::cout << "****** Forward kinematics example ******\n";
    std::cout << "****************************************\n";

    // Single input
    std::cout << "==========Single Input==========\n";
    double theta = M_PI * 130.0 / 180.0;
    double beta =  M_PI * 50.0 / 180.0;

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
    /* Contact map */
    legmodel.contact_map(theta, beta);
    std::cout << "Output rim with single value input: " << legmodel.rim << std::endl;
    std::cout << "Output alpha with single value input: " << legmodel.alpha << std::endl;

    /* Inverse kinematics */
    std::cout << "\n";
    std::cout << "****************************************" << std::endl;
    std::cout << "****** Inverse kinematics example ******" << std::endl;
    std::cout << "****************************************" << std::endl;
    // inverse for G
    std::cout << "==========Inverse for G==========" << std::endl;
    double G_p[2] = {0.05, -0.25};
    std::cout << "Input G: " << G_p[0] << ", " << G_p[1] << std::endl;
    legmodel.inverse(G_p, "G");
    std::cout << "Output theta, beta (degree): " << legmodel.theta*180.0/M_PI << ", "<< legmodel.beta*180.0/M_PI << std::endl;
    std::cout << "Output G: " << legmodel.G[0] << ", " << legmodel.G[1] << std::endl;
    // inverse for left upper rim
    std::cout << "==========Inverse for U_l==========" << std::endl;
    double Ul_p[2] = {-0.01, -0.015};
    std::cout << "Input U_l: " << Ul_p[0] << ", " << Ul_p[1] << std::endl;
    legmodel.inverse(Ul_p, "Ul");
    std::cout << "Output theta, beta (degree): " << legmodel.theta*180.0/M_PI << ", "<< legmodel.beta*180.0/M_PI << std::endl;
    std::cout << "Output U_l: " << legmodel.U_l[0] << ", " << legmodel.U_l[1] << std::endl;
    // inverse for right lower rim
    std::cout << "==========Inverse for L_r==========" << std::endl;
    double Lr_p[2] = {-0.01, -0.015};
    std::cout << "Input L_r: " << Lr_p[0] << ", " << Lr_p[1] << std::endl;
    legmodel.inverse(Lr_p, "Lr");
    std::cout << "Output theta, beta (degree): " << legmodel.theta*180.0/M_PI << ", "<< legmodel.beta*180.0/M_PI << std::endl;
    std::cout << "Output L_r: " << legmodel.L_r[0] << ", " << legmodel.L_r[1] << std::endl;



    // Eigen::VectorXd a(5);
    // Eigen::VectorXd b(5);
    // Eigen::VectorXd result(5);

    // // 初始化向量
    // a << 1.0, 2.0, 3.0, 4.0, 5.0;

    // b << 5.0, 4.0, 3.0, 2.0, 1.0;
    // a << 5.0, 2.0, 3.0, 4.0, 5.0;
    // std::cout << a ;


    
    auto start2 = std::chrono::high_resolution_clock::now();
    double angle;
    for (int i=0; i<10000000; i++){
        // angle = std::arg(std::complex<double>(1, 5) / std::complex<double>(0, -5));
        angle = std::atan2(5, 1) - std::atan2(-5, 0);
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    // 計算執行時間，並轉換成毫秒
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
    std::cout << "angle: " << angle  << ","<<std::atan2(5, 1) << ","<< std::atan2(-5, 0) << std::endl;
    std::cout << "time: " << duration2.count() << " ms" << std::endl;


    return 0;
}