#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <limits>
#include <cstddef>
#include <algorithm>


// 0. Defining operations and struct 

struct Vec3 {
    double x;
    double y;
    double z;
}; 

struct Grid2D {
    int nx;
    int ny; 

    std::vector<double> x;
    std::vector<double> y;

    std::vector<Vec3> magnetic_field;
    std::vector<Vec3> electric_field;
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

double dt = 0.01; 
double T_max = 300.0; 
int n_steps = static_cast<int>(std::ceil( T_max / dt)) + 1; 

int num_particles = 1000;
int random_seed = 42; 

Vec3 boris_push(Vec3 v, Vec3 B, Vec3 E, double q_over_m, double dt) {
    
    Vec3 v_minus = add(v, scalar(E, q_over_m * dt * 0.5)); 

    Vec3 t = scalar(B, q_over_m * dt * 0.5 ); 

    Vec3 v_prime = add(v_minus, cross(v_minus, t));

    Vec3 s = scalar(t, 2.0 / (1.0 + dot(t, t))); 

    Vec3 v_plus = add(v_minus, cross(v_prime, s));

    Vec3 v_next = add(v_plus, scalar(E, q_over_m * dt * 0.5)); 

    return v_next;
}

int main() {
    std::vector<Vec3> initial_velocities(num_particles);

    double mean_perp = 0.0049;
    double sigma_perp = 0.0049;
    double mean_par = 0.0105;
    double sigma_par = 0.0049;

    std::mt19937 generator(random_seed);
    std::normal_distribution<double> v_perp_distribution(
        mean_perp, sigma_perp
    );
    std::normal_distribution<double> v_par_distribution(
        mean_par, sigma_par
    );

    for (int p = 0; p < num_particles; ++p) {
        initial_velocities[p] = {
            v_perp_distribution(generator),
            v_perp_distribution(generator),
            v_par_distribution(generator)
        };
    }

    std::cout << "Number of particles: "
              << initial_velocities.size() << '\n';
    std::cout << "First velocity: "
              << initial_velocities[0].x << ", "
              << initial_velocities[0].y << ", "
              << initial_velocities[0].z << '\n';
}

