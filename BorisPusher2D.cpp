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
// 0.1. Defining the 3D Grid structure
struct Grid3D {
    int nx;
    int ny;
    int nz; 

    double x_min; 
    double y_min;
    double z_min;
    double spacing;

    std::vector<Vec3> magnetic_field;
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

Grid3D make_grid(double min_coord, double max_coord, double spacing, double coeff) {
    Grid3D grid{};
    grid.x_min = min_coord;
    grid.y_min = min_coord;
    grid.z_min = min_coord;
    grid.spacing = spacing;

    int points = static_cast<int>(std::round((max_coord - min_coord) / spacing)) + 1;

    grid.nx = points;
    grid.ny = points;
    grid.nz = points;

    std::size_t total_points = 
        static_cast<std::size_t>(points) * points * points;
    grid.magnetic_field.resize(total_points);

    for (int k = 0; k < points; ++k) {
    for (int j = 0; j < points; ++j) {
        for (int i = 0; i < points; ++i) {
            Vec3 position{
                min_coord + i * spacing,
                min_coord + j * spacing,
                min_coord + k * spacing
            };

            std::size_t index =
                (static_cast<std::size_t>(k) * points + j) * points + i;

            grid.magnetic_field[index] = dipole_B(position, coeff);
        }
    }
}
    return grid;
}

std::size_t grid_index(const Grid3D& grid, int i, int j, int k) {
    return (static_cast<std::size_t>(k) * grid.ny + j) * grid.nx + i;
}

Vec3 lerp(Vec3 a, Vec3 b, double fraction) {
    return add(
        scalar(a, 1.0 - fraction),
        scalar(b, fraction)
    );
}

bool interpolate_B(const Grid3D& grid, Vec3 position, Vec3& B) {
    double gx = (position.x - grid.x_min) / grid.spacing;
    double gy = (position.y - grid.y_min) / grid.spacing;
    double gz = (position.z - grid.z_min) / grid.spacing;

    if (gx < 0.0 || gy < 0.0 || gz < 0.0 ||
        gx > grid.nx - 1 || gy > grid.ny - 1 || gz > grid.nz - 1) {
            return false;
        }
    int i = static_cast<int>(std::floor(gx));
    int j = static_cast<int>(std::floor(gy));
    int k = static_cast<int>(std::floor(gz));

    double fx = gx - i;
    double fy = gy - j;
    double fz = gz - k;

    if (i == grid.nx - 1) { i --; fx = 1.0; }
    if (j == grid.ny - 1) { j --; fy = 1.0; }
    if (k == grid.nz - 1) { k --; fz = 1.0; }
    // Define each corner of the cube surrounding the point
    Vec3 c000 = grid.magnetic_field[grid_index(grid, i, j, k)];
    Vec3 c100 = grid.magnetic_field[grid_index(grid, i + 1, j, k)];
    Vec3 c010 = grid.magnetic_field[grid_index(grid, i, j + 1, k)];
    Vec3 c110 = grid.magnetic_field[grid_index(grid, i + 1, j + 1, k)];
    Vec3 c001 = grid.magnetic_field[grid_index(grid, i, j, k + 1)];
    Vec3 c101 = grid.magnetic_field[grid_index(grid, i + 1, j, k + 1)];
    Vec3 c011 = grid.magnetic_field[grid_index(grid, i, j + 1, k + 1)];
    Vec3 c111 = grid.magnetic_field[grid_index(grid, i + 1, j + 1, k + 1)];
    // Interpolate along 
    Vec3 c00 = lerp(c000, c100, fx);
    Vec3 c10 = lerp(c010, c110, fx);
    Vec3 c01 = lerp(c001, c101, fx);
    Vec3 c11 = lerp(c011, c111, fx);
    // interpolate along y 
    Vec3 c0 = lerp(c00, c10, fy);
    Vec3 c1 = lerp(c01, c11, fy);
    // Interpolate along z
    B = lerp(c0, c1, fz);
    return true;
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
// 3. Create a grid and sanity check
Grid3D dipole_grid = make_grid(-2.5, 2.5, 0.025, dipole_coeff);
std::cout << "Grid points per axis: " << dipole_grid.nx << '\n'; 
std::cout << "Stored B values: " << dipole_grid.magnetic_field.size() << '\n';
// 4. 
Vec3 test_position{1.237, -0.483, 0.918};
Vec3 B_from_grid; 

bool inside_grid = interpolate_B(dipole_grid, test_position, B_from_grid);

if (!inside_grid) {
    std::cout << "Test position is outside the grid.\n";
    return 1;
}

Vec3 B_exact = dipole_B(test_position, dipole_coeff);
Vec3 field_error = subtract(B_from_grid, B_exact);

double relative_error = std::sqrt(dot(field_error, field_error)) / std::sqrt(dot(B_exact, B_exact));

std::cout << "Analytic B: "
          << B_exact.x << ", " << B_exact.y << ", " << B_exact.z << '\n';
std::cout << "Grid B:     "
          << B_from_grid.x << ", " << B_from_grid.y << ", " << B_from_grid.z << '\n';
std::cout << "Relative field error: " << relative_error << '\n';

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
std::vector<unsigned char> grid_is_lost(num_particles, false);
std::vector<double> grid_hit_times(num_particles, std::numeric_limits<double>::quiet_NaN());
std::vector<unsigned char> grid_out_of_bounds(num_particles, false);

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

#pragma omp parallel for   
for (int p = 0; p < num_particles; ++p) {
    Vec3 r = initial_positions[p];
    Vec3 v_at_zero = initial_velocities[p]; 

    if (dot(r, r) < r_earth * r_earth) {
        grid_is_lost[p] = true;
        grid_hit_times[p] = 0.0;
        continue;
    }

    Vec3 B0; 
    if (!interpolate_B(dipole_grid, r, B0)) {
        grid_out_of_bounds[p] = true;
        continue;
    }

    Vec3 v_half = boris_rotate(v_at_zero, B0, q_over_m, -0.5 * dt);

    for (int i = 1; i < n_steps; ++i) {
        Vec3 B;
        if (!interpolate_B(dipole_grid, r, B)) {
            grid_out_of_bounds[p] = true;
            break;
        }

        Vec3 v_next_half = boris_rotate(v_half, B, q_over_m, dt);
        Vec3 r_next = add(r, scalar(v_next_half, dt));

        double fraction = segment_sphere_fraction(r, r_next, r_earth);

        if (fraction >= 0.0) {
            grid_is_lost[p] = true;
            grid_hit_times[p] = (i - 1 + fraction) * dt;
            break;
        }

        r = r_next;
        v_half = v_next_half;
    }
}

int function_impact_count = 0;
int grid_impact_count = 0;
int grid_outside_count = 0;
int outcome_mismatches = 0;

std::vector<double> hit_time_differences;

for (int p = 0; p < num_particles; ++p) {
    if (is_lost[p]) {
        ++function_impact_count;
    }

    if (grid_is_lost[p]) {
        ++grid_impact_count;
    }

    if (grid_out_of_bounds[p]) {
        ++grid_outside_count;
        continue;
    }

    if (is_lost[p] != grid_is_lost[p]) {
        ++outcome_mismatches;
    }

    if (is_lost[p] && grid_is_lost[p]) {
        double difference = std::abs(hit_times[p] - grid_hit_times[p]);
        hit_time_differences.push_back(difference);
    }
}

std::cout << "Analytic dipole impacts: " << function_impact_count << '\n';
std::cout << "Grid dipole impacts:     " << grid_impact_count << '\n';
std::cout << "Particles leaving grid:  " << grid_outside_count << '\n';
std::cout << "Impact mismatches:       " << outcome_mismatches << '\n';

if (!hit_time_differences.empty()) {
    std::sort(hit_time_differences.begin(), hit_time_differences.end());

    std::size_t n = hit_time_differences.size();
    double mean_difference = 0.0;

    for (double difference : hit_time_differences) {
        mean_difference += difference;
    }
    mean_difference /= static_cast<double>(n);

    double median_difference;
    if (n % 2 == 0) {
        median_difference =
            (hit_time_differences[n / 2 - 1] +
             hit_time_differences[n / 2]) / 2.0;
    } else {
        median_difference = hit_time_differences[n / 2];
    }

    std::size_t p95_index =
        static_cast<std::size_t>(std::ceil(0.95 * n)) - 1;

    std::cout << "Impacts in both runs: " << n << '\n';
    std::cout << "Mean hit-time difference: " << mean_difference << '\n';
    std::cout << "Median hit-time difference: " << median_difference << '\n';
    std::cout << "95th percentile difference: "
              << hit_time_differences[p95_index] << '\n';
    std::cout << "Maximum hit-time difference: "
              << hit_time_differences.back() << '\n';
}
}
