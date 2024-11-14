#include <iostream>
#include <vector>
#include <cmath>

// Polynomial class to represent and evaluate polynomials
class Polynomial {
public:
    std::vector<double> coefficients;

    Polynomial(const std::vector<double>& coefs) : coefficients(coefs) {}

    // Evaluate polynomial at x
    double evaluate(double x) const {
        double result = 0.0;
        double power = 1.0;
        for (double coef : coefficients) {
            result += coef * power;
            power *= x;
        }
        return result;
    }

    // Derivative of the polynomial
    Polynomial derivative() const {
        std::vector<double> deriv_coefs;
        for (size_t i = 1; i < coefficients.size(); ++i) {
            deriv_coefs.push_back(coefficients[i] * i);
        }
        return Polynomial(deriv_coefs);
    }
};

// Define polynomial coefficients
std::vector<double> A_x_coef = {1.1298507215688325e-05, -0.08009406847959913, 0.00031078052601928794, 0.012797805183985959, 0.0005296240879456367, -0.0009747108437109868, 0.00010114657513265072, 3.9686108833838617e-07};
std::vector<double> A_y_coef = {0.07999814821731825, 1.7151920935796158e-05, -0.040064720911549265, 0.00013130348527951678, 0.003174527009668425, 0.00011925178438898585, -0.00016662716127708877, 1.5155513490796227e-05};

// Create polynomial objects
Polynomial A_x_poly(A_x_coef);
Polynomial A_y_poly(A_y_coef);

// Example usage
int main() {
    double x = 0.5; // Example input
    std::cout << "A_x evaluated at " << x << ": " << A_x_poly.evaluate(x) << std::endl;
    std::cout << "A_y evaluated at " << x << ": " << A_y_poly.evaluate(x) << std::endl;

    // Compute derivatives
    Polynomial A_x_deriv = A_x_poly.derivative();
    Polynomial A_y_deriv = A_y_poly.derivative();
    std::cout << "A_x derivative evaluated at " << x << ": " << A_x_deriv.evaluate(x) << std::endl;
    std::cout << "A_y derivative evaluated at " << x << ": " << A_y_deriv.evaluate(x) << std::endl;
    std::cout << A_y_deriv.coefficients[0];
    return 0;
}
