
/**
 * This basic Kalman Filter demo simulate a moving object with constant
 * velocity and random acceleration noise.
 *
 * The system is described with equation x(k) = A * x(k-1) + B * u(k-1)
 * The output is described with equation z(k) = H * x(k)
 */
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "kalman_filter.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// Time and samples
constexpr double SYSTEM_DT = 0.01;     // The system works with rate of 100 Hz
constexpr double MEASUREMENT_DT = 0.1; // The measurements come with rate of 10 Hz
constexpr double SIMULATION_TIME = 10; // The time of simulation is 10 s

constexpr int N = static_cast<int>(SIMULATION_TIME / SYSTEM_DT); // simulation duration
constexpr int M = static_cast<int>(MEASUREMENT_DT / SYSTEM_DT);  // measurement duration

int main([[maybe_unused]]int argc,
         [[maybe_unused]]char *argv[]) {

    // Buffers for plots
    std::vector<double> time(N);

    std::vector<double> zero_pos(N, 0);
    std::vector<double> true_pos_u(N);
    std::vector<double> true_pos_v(N);
    std::vector<double> true_vel_u(N);
    std::vector<double> true_vel_v(N);

    std::vector<double> measured_pos_u(N);
    std::vector<double> measured_pos_v(N);
    std::vector<double> estimated_pos_u(N);
    std::vector<double> estimated_pos_v(N);

    // Pseudo random numbers generator
    double measurement_mu = 0.0;     // Mean
    double measurement_sigma = 0.05; // Standard deviation

    double process_mu = 0.0;    // Mean
    double process_sigma = 2.0; // Standard deviation

    std::default_random_engine generator;
    std::normal_distribution<double> measurement_noise(measurement_mu, measurement_sigma);
    std::normal_distribution<double> process_noise(process_mu, process_sigma);
    // clang-format off
    // Preparing KF
    // x = [Pu, Pv, Vu, Vv]^T
    Eigen::Matrix4d A; /* 4x4 */
    A << 1.0, 0.0, SYSTEM_DT, 0.0,
         0.0, 1.0, 0.0, SYSTEM_DT,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;

    Eigen::Vector4d B;
    B.setZero();

    Eigen::Matrix<double, 2, 4> H; /* 2x4 */
    H << 1.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0;

    Eigen::Matrix4d P0;
    P0.setIdentity();

    // The process and measurement covariances are sort of tunning parameters
    Eigen::Matrix4d Q; /* 4x4 */
    double sigma_au_2 = 0.64;
    double sigma_av_2 = 0.64;
    double dt_2 = pow(SYSTEM_DT, 2);
    double dt_3 = pow(SYSTEM_DT, 3);
    double dt_4 = pow(SYSTEM_DT, 4);

    double q11 = dt_4 * sigma_au_2 / 4;
    double q13 = dt_3 * sigma_au_2 / 2;
    double q22 = dt_4 * sigma_av_2 / 4;
    double q24 = dt_3 * sigma_av_2 / 2;
    double q31 = dt_3 * sigma_au_2 / 2;
    double q33 = dt_2 * sigma_au_2;
    double q42 = dt_3 * sigma_av_2 / 2;
    double q44 = dt_2 * sigma_av_2;
    Q << q11, 0,   q13, 0,
         0,   q22, 0,   q24,
         q31, 0,   q33, 0,
         0,   q42, 0,   q44;

    Eigen::Matrix2d R; /* 2x2 */
    R << 0.01, 0.0,
         0.0,  0.01;

    Eigen::Matrix<double, 4, 2> G; /* 4x2 */
    G << 1, 0,
         0, 1,
         0, 0,
         0, 0;

    // clang-format on
    kf::KalmanFilter filter(A, B, H, Q, R, P0);
    filter.init();

    // Initial values (unknown by KF)
    time[0] = 0.0;
    true_pos_u[0] = 1.0;
    true_vel_u[0] = 0.1;
    measured_pos_u[0] = 0.0;
    estimated_pos_u[0] = 0.0;

    true_pos_v[0] = 0.0;
    true_vel_v[0] = 0.5;
    measured_pos_v[0] = 0.0;
    estimated_pos_v[0] = 0.0;

    Eigen::VectorXd z(2); /* 2x1 */
    Eigen::VectorXd e(2); /* 2x1 */
    double l_threshold = 20.0;

    // Simulation
    for (size_t i = 1; i < N; ++i) {
        time[i] = i * SYSTEM_DT;

        if (i == N / 2) {
            // switch target and set new target state
            true_vel_u[i] = 0.3;
            true_pos_u[i] = 3.0;

            true_vel_v[i] = 0.3;
            true_pos_v[i] = 3.0;
        } else {
            true_vel_u[i] = true_vel_u[i - 1] + process_noise(generator) * SYSTEM_DT;
            true_pos_u[i] = true_pos_u[i - 1] + true_vel_u[i] * SYSTEM_DT +
                            0.5 * process_noise(generator) * dt_2;

            true_vel_v[i] = true_vel_v[i - 1] + process_noise(generator) * SYSTEM_DT;
            true_pos_v[i] = true_pos_v[i - 1] + true_vel_v[i] * SYSTEM_DT +
                            0.5 * process_noise(generator) * dt_2;
        }

        // New measurement comes once every M samples of the system
        if (i % M == 1) {
            measured_pos_u[i] = true_pos_u[i] + measurement_noise(generator);
            measured_pos_v[i] = true_pos_v[i] + measurement_noise(generator);
        } else {
            measured_pos_u[i] = measured_pos_u[i - 1];
            measured_pos_v[i] = measured_pos_v[i - 1];
        }

        z(0) = measured_pos_u[i];
        z(1) = measured_pos_v[i];

        filter.predict();

        /* chi-square test */
        e = z - H * filter.getPredict();
        double r = e.transpose() * filter.getInnovationCovariance().inverse() * e;
        if (r > l_threshold) {
            printf("[Kalman Filter]: Switch target! threshold: %.3f, real: %.3f\n", l_threshold, r);
            filter.setEstimateCovariance(P0);
            filter.setEstimate(G * z);
            estimated_pos_u[i] = z(0);
            estimated_pos_v[i] = z(1);
            continue;
        }

        filter.update(z);

        estimated_pos_u[i] = filter.getEstimate()(0);
        estimated_pos_v[i] = filter.getEstimate()(1);
    }

    // Plot
#if KF_USE_MATPLOTLIB
    plt::figure_size(1280, 768);
    plt::title("Estimate of position");
    plt::xlabel("Time [s]");
    plt::ylabel("Position [m]");
    plt::named_plot("Truth: u", time, true_pos_u, "--");
    plt::named_plot("Measure: u", time, measured_pos_u, "-");
    plt::named_plot("Estimate: u", time, estimated_pos_u, "-");
    plt::named_plot("Truth: v", time, true_pos_v, "--");
    plt::named_plot("Measure: v", time, measured_pos_v, "-");
    plt::named_plot("Estimate: v", time, estimated_pos_v, "-");
    plt::legend();
    plt::grid(true, {{"linestyle", "--"}});
    plt::xlim(0.0, SIMULATION_TIME - SYSTEM_DT);
    plt::save("assets/kalman_filter_acc_noise_model_chi-square.png");
    plt::show();
#endif

    return 0;
}