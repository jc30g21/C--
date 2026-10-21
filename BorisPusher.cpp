#include <iostream>
#include <cmath>

// 0. Defining operations and struct 

struct Vec3 {
    double x;
    double y;
    double z;
}; 

Vec3 add(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z}; 
}

double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(Vec3 a, Vec3 b) { 
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

Vec3 scalar(Vec3 a, double b) {
    return {a.x * b, a.y * b, a.z * b};
}

constexpr double pi = 3.14159265358979323846;

// 1. Natural units and parameters

double q = 1.0; 
double m = 1.0; 
double q_over_m = q / m; 

double mu0 = 1.0; 
double dipole_moment = 300.0; 
double dipole_coeff = mu0 * dipole_moment / (4.0 * pi); 

double r_earth = 1.0; 
double l_shell = 2.0;

double dt = 0.01; 
double T_max = 300.0; 
int n_steps = static_cast<int>(std::ceil( T_max / dt)) + 1; 

int num_particles = 1000;
int random_seed = 42; 

// 2. 

Vec3 initial_position{l_shell, 0, 0};





int main() {
    Vec3 position{2.0, 0.0, 0.0};
    Vec3 displacement{0.1, 0.2, 0.3};

    Vec3 new_position = add(position, displacement);

    Vec3 difference = subtract(position, displacement);

    std::cout << "New position: "
              << new_position.x << ", "
              << new_position.y << ", "
              << new_position.z << '\n'; 

    std::cout << "Difference: "
              << difference.x << ", "
              << difference.y << ", "
              << difference.z << '\n';

    std::cout << "Number of steps: " << n_steps << '\n';

    std::cout << "Initial position: "
          << initial_position.x << ", "
          << initial_position.y << ", "
          << initial_position.z << '\n';
}