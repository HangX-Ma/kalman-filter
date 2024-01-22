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
#include <fmt/core.h>

namespace kf {

class KalmanFilter {
    public:
        /**
         * Create a Kalman filter with the specified matrices.
         *   A - System dynamics matrix
         *   H - Measurement matrix
         *   Q - Process noise covariance
         *   R - Measurement noise covariance
         *   P - Estimate error covariance
         */
      KalmanFilter(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B, const Eigen::MatrixXd &H,
                   const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R, const Eigen::MatrixXd &P)
            : A_(A), B_(B), H_(H),
              l_(B_.cols()), m_(H_.rows()), n_(A_.rows()),
              P_(P), P0_(P), Q_(Q), R_(R)
        {
            K_.resize(n_, m_);
            I_.resize(n_, n_);

            P_.resize(n_, n_);
            S_.resize(m_, m_);

            x_pred_.resize(n_);
            x_est_.resize(n_);
            u_.resize(l_);

            K_.setZero();
            I_.setIdentity();

            S_.setIdentity();
        }

        void init(const Eigen::VectorXd &x0) {
            if (x0.size() != x_est_.size()) {
                throw std::length_error(fmt::format("Incorrect dimensions of estimated state vector in {}", __func__));
            }
            x_est_ = x0;
            P_ = P0_;
        }

        void init() {
            x_est_.setZero();
            P_ = P0_;
        }

        void setEstimate(const Eigen::VectorXd & x_est) {
            if (x_est.size() != x_est_.size()) {
                throw std::length_error(fmt::format("Incorrect dimensions of estimated state vector in {}", __func__));
            }
            x_est_ = x_est;
        }

        void setInput(const Eigen::VectorXd& u) {
            if (u.size() != u_.size()) {
                throw std::length_error(fmt::format("Incorrect dimensions of input vector in {}", __func__));
            }
            u_ = u;
        }

        void setEstimateCovariance(const Eigen::MatrixXd& P) {
            if (P.size() != P_.size()) {
                throw std::length_error(fmt::format("Incorrect dimensions of estimate cov. matrix in {}", __func__));
            }
            P_ = P;
        }

        void setProcessCovariance(const Eigen::MatrixXd& Q) {
            if (Q.size() != Q_.size()) {
                throw std::length_error(fmt::format("Incorrect dimensions of process cov. matrix in {}", __func__));
            }
            Q_ = Q;
        }

        /** @brief Estimate covariance matrix getter
         *
         * @returns the current estimate covariance matrix P
         */
        const Eigen::MatrixXd& getEstimateCovariance(void) const { return P_; }


        /** @brief Innovation covariance matrix getter
         *
         * Innovation covariance matrix is calculated during correction step
         * S = H * P * trans(H) + R.
         *
         * @returns the current innovation covariance matrix S
         */
        const Eigen::MatrixXd& getInnovationCovariance(void) const { return S_; }

        const Eigen::VectorXd& getPredict(void) const { return x_pred_; };
        const Eigen::VectorXd& getEstimate(void) const { return x_est_; };

        /** @brief Performs the KF prediction step
         *
         * Calculates the state predicted from the model and updates the estimate
         * error covariance matrix.
         */
        void predict(void) {
            x_pred_ = A_ * x_est_ + B_ * u_;
            P_ = A_ * P_ * A_.transpose() + Q_;
            // prepare for 'update'
            S_ = H_ * P_ * H_.transpose() + R_;
        }

        void predict(const Eigen::VectorXd& u) {
            setInput(u);
            predict();
        }

        /** @brief Performs the KF correction step
         *
         * Calculates the new Kalman gain, corrects the state prediction to obtain new
         * state estimate, and updates estimate error covariance as well as innovation
         * covariance.
         */
        void update(const Eigen::VectorXd z) {
            K_ = P_ * H_.transpose() * S_.inverse();
            x_est_ = x_pred_ + K_ * (z - H_ * x_pred_);
            P_ = (I_ - K_ * H_) * P_;
        }

    private:
        // Matrices for computation
        Eigen::MatrixXd A_, B_, H_;

        size_t l_; // dimension of input vector. Can be zero for autonomous system.
        size_t m_; // dimension of output vector. Can be zero for prediction only.
        size_t n_; // dimension of state vector

        // n-size identity
        Eigen::MatrixXd I_;

        // Kalman gain matrix:
        Eigen::MatrixXd K_;     /* n x m */

        // Covariance matrices:
        Eigen::MatrixXd P_; // estimate, n x n
        Eigen::MatrixXd P0_; // initial estimate value, n x n
        Eigen::MatrixXd Q_; // process noise, n x n
        Eigen::MatrixXd R_; // measurement noise, m x m
        Eigen::MatrixXd S_; // innovation, m x m

        // input vector
        Eigen::VectorXd u_;

        // Estimated states
        Eigen::VectorXd x_pred_, x_est_;
};

}

#endif  //!__KALMAN_FILTER__H__