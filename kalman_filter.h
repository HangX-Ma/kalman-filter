/**
 * @file kalman_filter.h
 * @author HangX-Ma (contour.9x@gmail.com)
 * @brief basic kalman filter implementation
 * @version 0.1
 * @date 2023-08-30
 */

#ifndef __KALMAN_FILTER__H__
#define __KALMAN_FILTER__H__

#include <Eigen/Eigen>
#include <utility>
#include <fmt/core.h>

namespace kf
{

class KalmanFilter
{
public:
    /**
     * Create a Kalman filter with the specified matrices.
     *   A - System dynamics matrix
     *   H - Measurement matrix
     *   Q - Process noise covariance
     *   R - Measurement noise covariance
     *   P - Estimate error covariance
     */
    KalmanFilter(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd H, Eigen::MatrixXd Q,
                 Eigen::MatrixXd R, const Eigen::MatrixXd &P)
        : A(std::move(A)), B(std::move(B)), H(std::move(H)), l(B.cols()), m(H.rows()), n(A.rows()),
          P(P), P0(P), Q(std::move(Q)), R(std::move(R))
    {
        K.resize(n, m);
        I.resize(n, n);
        S.resize(m, m);

        x_pred.resize(n);
        x_est.resize(n);
        u.resize(l);

        K.setZero();
        I.setIdentity();
        S.setIdentity();
    }

    void init(const Eigen::VectorXd &x0)
    {
        if (x0.size() != x_est.size()) {
            throw std::length_error(
                fmt::format("Incorrect dimensions of estimated state vector in {}", __func__));
        }
        x_est = x0;
        P = P0;
    }

    void init()
    {
        x_est.setZero();
        P = P0;
    }

    void setEstimate(const Eigen::VectorXd &x_est)
    {
        if (x_est.size() != this->x_est.size()) {
            throw std::length_error(
                fmt::format("Incorrect dimensions of estimated state vector in {}", __func__));
        }
        this->x_est = x_est;
    }

    void setInput(const Eigen::VectorXd &u)
    {
        if (u.size() != this->u.size()) {
            throw std::length_error(
                fmt::format("Incorrect dimensions of input vector in {}", __func__));
        }
        this->u = u;
    }

    void setEstimateCovariance(const Eigen::MatrixXd &P)
    {
        if (P.size() != this->P.size()) {
            throw std::length_error(
                fmt::format("Incorrect dimensions of estimate cov. matrix in {}", __func__));
        }
        this->P = P;
    }

    void setProcessCovariance(const Eigen::MatrixXd &Q)
    {
        if (Q.size() != this->Q.size()) {
            throw std::length_error(
                fmt::format("Incorrect dimensions of process cov. matrix in {}", __func__));
        }
        this->Q = Q;
    }

    /** @brief Estimate covariance matrix getter
     *
     * @returns the current estimate covariance matrix P
     */
    const Eigen::MatrixXd &getEstimateCovariance() const { return P; }

    /** @brief Innovation covariance matrix getter
     *
     * Innovation covariance matrix is calculated during correction step
     * S = H * P * trans(H) + R.
     *
     * @returns the current innovation covariance matrix S
     */
    const Eigen::MatrixXd &getInnovationCovariance() const { return S; }

    const Eigen::VectorXd &getPredict() const { return x_pred; };
    const Eigen::VectorXd &getEstimate() const { return x_est; };

    /** @brief Performs the KF prediction step
     *
     * Calculates the state predicted from the model and updates the estimate
     * error covariance matrix.
     */
    void predict()
    {
        x_pred = A * x_est + B * u;
        P = A * P * A.transpose() + Q;
        // prepare for 'update'
        S = H * P * H.transpose() + R;
    }

    void predict(const Eigen::VectorXd &u)
    {
        setInput(u);
        predict();
    }

    /** @brief Performs the KF correction step
     *
     * Calculates the new Kalman gain, corrects the state prediction to obtain new
     * state estimate, and updates estimate error covariance as well as innovation
     * covariance.
     */
    void update(const Eigen::VectorXd &z)
    {
        K = P * H.transpose() * S.inverse();
        x_est = x_pred + K * (z - H * x_pred);
        P = (I - K * H) * P;
    }

protected:
    // Matrices for computation
    Eigen::MatrixXd A, B, H;

    size_t l; // dimension of input vector. Can be zero for autonomous system.
    size_t m; // dimension of output vector. Can be zero for prediction only.
    size_t n; // dimension of state vector

    // n-size identity
    Eigen::MatrixXd I;

    // Kalman gain matrix:
    Eigen::MatrixXd K; /* n x m */

    // Covariance matrices:
    Eigen::MatrixXd P;  // estimate, n x n
    Eigen::MatrixXd P0; // initial estimate value, n x n
    Eigen::MatrixXd Q;  // process noise, n x n
    Eigen::MatrixXd R;  // measurement noise, m x m
    Eigen::MatrixXd S;  // innovation, m x m

    // input vector
    Eigen::VectorXd u;

    // Estimated states
    Eigen::VectorXd x_pred, x_est;
};

} // namespace kf

#endif //!__KALMAN_FILTER__H__