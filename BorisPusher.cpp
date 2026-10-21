#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <limits>

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
double dipole_coeff = mu0 * sdipole_moment / (4.0 * pi); 

double r_earth = 1.0; 
double l_shell = 2.0;

double dt = 0.01; 
double T_max = 300.0; 
int n_steps = static_cast<int>(std::ceil( T_max / dt)) + 1; 

int num_particles = 1000000;
int random_seed = 42; 

Vec3 dipole_B(Vec3 position, double coeff) {
    double r2 = dot(position, position);

    if (r2 < 1.0e-24) {
        return {0.0, 0.0, 0.0};
    }
    double r = std::sqrt(r2);
    double inv_r5 = 1.0 / (r2 * r2 * r);

    return {
        coeff * (3.0 * position.z * position.x * inv_r5),
        coeff * (3.0 * position.z * position.y * inv_r5),
        coeff * ((3.0 * position.z * position.z - r2) * inv_r5)
    };
}

Vec3 boris_rotate(Vec3 v, Vec3 B, double q_over_m, double dt) {
    Vec3 t = scalar(B, q_over_m * dt * 0.5); 

    Vec3 v_prime = add(v, cross(v, t));

    double t_squared = dot(t, t);

    Vec3 s = scalar(t, 2.0 / (1.0 + t_squared));

    Vec3 v_next = add(v, cross(v_prime, s)); 

    return v_next;
}

double segment_sphere_fraction(Vec3 r0, Vec3 r1, double radius) {
    Vec3 delta = subtract(r1, r0);
    double a = dot(delta, delta);

    if (a == 0.0) {
        return -1.0;
    }

    double b = 2.0 * dot(r0, delta);
    double c = dot(r0, r0) - radius * radius;
    double discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) {
        return -1.0;
    }

    double root = std::sqrt(discriminant);
    double fraction1 = (-b - root) / (2.0 * a);
    double fraction2 = (-b + root) / (2.0 * a);

    if (0.0 <= fraction1 && fraction1 <= 1.0) {
        return fraction1;
    }

    if (0.0 <= fraction2 && fraction2 <= 1.0) {
        return fraction2;
    }

    return -1.0;
}

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

std::vector<unsigned char> is_lost(num_particles, false);
std::vector<double> hit_times(num_particles, std::numeric_limits<double>::quiet_NaN());

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

Vec3 B_at_equator = dipole_B(Vec3{l_shell, 0.0, 0.0}, dipole_coeff);
std::cout << "Magnetic field at equator: "
          << B_at_equator.x << ", "
          << B_at_equator.y << ", "
          << B_at_equator.z << '\n';

Vec3 test_v{0.2, 0.3, 0.4};
Vec3 test_B{0.0, 0.0, 2.0};

Vec3 rotated_v = boris_rotate(test_v, test_B, q_over_m, dt);

double speed_before = std::sqrt(dot(test_v, test_v));
double speed_after = std::sqrt(dot(rotated_v, rotated_v));

std::cout << "Test speed before: " << speed_before << '\n';
std::cout << "Test speed after:  " << speed_after << '\n';

#pragma omp parallel for
for (int p = 0; p < num_particles; ++p) {
    Vec3 r = initial_positions[p];
    Vec3 v_at_zero = initial_velocities[p];
    Vec3 v_half = boris_rotate(
        v_at_zero, dipole_B(r, dipole_coeff), q_over_m, -0.5 *dt
    ); 

    if (dot(r, r) < r_earth * r_earth) {
        is_lost[p] = true;
        hit_times[p] = 0.0;
        continue;
    }

    for (int i = 1; i < n_steps; ++i) { 
        Vec3 B = dipole_B(r, dipole_coeff);
        Vec3 v_next_half = boris_rotate(v_half, B, q_over_m, dt);
        Vec3 r_next = add(r, scalar(v_next_half, dt));

        double fraction = segment_sphere_fraction(r, r_next, r_earth);

        if (fraction >= 0.0) {
            is_lost[p] = true;
            hit_times[p] = ( i - 1 + fraction) * dt;
            break;
        }

        r = r_next;
        v_half = v_next_half;
    }
}
int impact_count = 0;

for (int p = 0; p < num_particles; ++p) {
    if (is_lost[p]) {
        ++impact_count;
    }
}

std::cout << "Impacts by T=" << T_max << ": "
          << impact_count << '\n';
}
