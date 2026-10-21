#include <iostream>
#include <cmath>
#include <vector>
#include <random>


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

// 2. Initial velocity and VDF

int main() {
    std::vector<Vec3> initial_positions(num_particles, Vec3{l_shell, 0.0, 0.0});
    std::cout << "Number of particles: " << initial_positions.size() << '\n';
std::vector<Vec3> initial_velocities(num_particles);
double mean_perp = 0.0049;
double sigma_perp = 0.0049;
double mean_par = 0.0105;
double sigma_par = 0.0049;

std::mt19937 generator(random_seed);
std::normal_distribution<double> v_per1_distribution(mean_perp, sigma_perp);
std::normal_distribution<double> v_par_distribution(mean_par, sigma_par);

for (int p = 0; p < num_particles; ++p) {
    initial_velocities[p] = {
        v_per1_distribution(generator),
        v_per1_distribution(generator),
        v_par_distribution(generator)
    };
}

std::vector<double> initial_v_perp(num_particles);
std::vector<double> initial_v_par(num_particles);

for (int p = 0; p < num_particles; ++p) {
    const Vec3& v = initial_velocities[p];
    initial_v_perp[p] = std::sqrt(v.x * v.x + v.y * v.y);
    initial_v_par[p] = std::abs(v.z);
}

std::cout << "First velocity: "
          << initial_velocities[0].x << ", "
          << initial_velocities[0].y << ", "
          << initial_velocities[0].z << '\n';  

double B_eq = dipole_coeff / (l_shell * l_shell * l_shell);
double R_curvature_eq = l_shell / 3.0;
double omega_eq = std::abs(q_over_m * B_eq);
double B_max_at_surface = 2.0 * dipole_coeff / (r_earth * r_earth * r_earth);
double omega_max_at_surface = std::abs(q_over_m * B_max_at_surface);

std::vector<double> rho_over_Rc(num_particles);

for (int p = 0; p < num_particles; ++p) {
    rho_over_Rc[p] = initial_v_perp[p] / (R_curvature_eq * omega_eq);}

std::cout << "Equatorial field at L=" << l_shell << ": " << B_eq << '\n';
std::cout << "Injection rho/Rc: first particle = " << rho_over_Rc[0] << '\n';
std::cout << "Maximum omega*dt at surface: " << omega_max_at_surface * dt << '\n';

if (omega_max_at_surface * dt > 0.05) {
    std::cout << "WARNING: omega*dt at the surface exceeds 0.05. "
              << "Consider reducing dt.\n";
}
}
