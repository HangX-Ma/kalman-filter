/**
 * This basic Kalman Filter demo simulate a moving object with constant
 * velocity.
 *
 * The system is described with equation x(k) = A * x(k-1) + B * u(k-1)
 * The output is described with equation z(k) = H * x(k)
*/
#include <string>
#include <vector>
#include <random>

#include "kalman_filter.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// Time and samples
const double system_dt = 0.01;        // The system works with rate of 100 Hz
const double measurement_dt = 0.1;    // The measurements come with rate of 10 Hz
const double simulation_time = 10;    // The time of simulation is 10 s

const int N = static_cast<int>(simulation_time / system_dt);    // simulation duration
const int M = static_cast<int>(measurement_dt / system_dt);     // measurement duration

int main (int argc,char *argv[]) {
    (void)argc;
    (void)argv;

    // Buffers for plots
    std::vector<double> time(N);

    std::vector<double> zero_pos(N, 0);
    std::vector<double> true_pos(N);
    std::vector<double> true_vel(N);

    std::vector<double> measured_pos(N);
    std::vector<double> estimated_pos(N);

    // Pseudo random numbers generator
    double measurement_mu = 0.0;      // Mean
    double measurement_sigma = 0.02;  // Standard deviation

    double process_mu = 0.0;          // Mean
    double process_sigma = 0.005;      // Standard deviation

    std::default_random_engine generator;
    std::normal_distribution<double> measurement_noise(measurement_mu, measurement_sigma);
    std::normal_distribution<double> process_noise(process_mu, process_sigma);

    // Preparing KF
    Eigen::Matrix2d A; /* 2x2 */
    A << 1.0, system_dt,
         0.0, 1.0;

    Eigen::Vector2d B;
    B.setZero();

    Eigen::RowVector2d H; /* 1x2 */
    H << 1.0, 0.0;

    // The process and measurement covariances are sort of tunning parameters
    Eigen::Matrix2d Q; /* 2x2 */
    Q << 0.001, 0.0,
         0.0,   0.001;

    Eigen::MatrixXd R(1, 1); /* 1x1 */
    R << 1.0;

    kf::KalmanFilter filter(A, B, H, Q, R);

    // Initial values (unknown by KF)
    time[0] = 0.0;
    true_pos[0] = 1.0;
    true_vel[0] = 0.1;

    measured_pos[0] = 0.0;
    estimated_pos[0] = 0.0;

    Eigen::VectorXd z(1); /* 1x1 */
    // Simulation
    for (size_t i = 1; i < N; ++i) {
        time[i] = i * system_dt;

        true_vel[i] = true_vel[i - 1] + process_noise(generator);
        true_pos[i] = true_pos[i - 1] + true_vel[i] * system_dt;

        // New measurement comes once every M samples of the system
        if (i % M == 1) {
            measured_pos[i] = true_pos[i] + measurement_noise(generator);
        } else {
            measured_pos[i] = measured_pos[i - 1];
        }
        z(0) = measured_pos[i];

        filter.predict();
        filter.update(z);

        estimated_pos[i] = filter.getEstimate()(0);
    }

    // Plot
#if KF_USE_MATPLOTLIB
    plt::figure_size(1280, 768);
    plt::title("Estimate of position");
    plt::xlabel("Time [s]");
    plt::ylabel("Position [m]");
    plt::named_plot("Truth", time, true_pos, "--");
    plt::named_plot("Measure", time, measured_pos, "-");
    plt::named_plot("Estimate", time, estimated_pos, "-");
    plt::legend();
    plt::grid(true, {{"linestyle", "--"}});
    plt::xlim(0.0, simulation_time - system_dt);
    plt::save("assets/kalman_filter_const_model.png");
    plt::show();
#endif

    return 0;
}