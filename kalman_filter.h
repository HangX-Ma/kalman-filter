/**
 * @file kalman_filter.h
 * @author HangX-Ma (contour.9x@gmail.com)
 * @brief basic kalman filter implementation
 * @version 0.1
 * @date 2023-08-30
 */

#ifndef __KALMAN_FILTER__H__
#define __KALMAN_FILTER__H__

#include "Eigen/Eigen"

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
        KalmanFilter(const Eigen::MatrixXd& A,
                     const Eigen::MatrixXd& B,
                     const Eigen::MatrixXd& H,
                     const Eigen::MatrixXd& Q,
                     const Eigen::MatrixXd& R)
            : A_(A), B_(B), H_(H), Q_(Q), R_(R),
              l_(B_.cols()), m_(H_.rows()), n_(A_.rows())
        {
            I_.resize(n_, n_);
            P_.resize(n_, n_);
            K_.resize(n_, m_);
            x_pred_.resize(n_);
            x_est_.resize(n_);
            u_.resize(l_);

            I_.setIdentity();
            P_.setIdentity();
            K_.setZero();
        }

        void setInput(const Eigen::VectorXd& u) {
            if (u.size() != u_.size()) {
                throw std::length_error("Incorrect dimensions of input vector");
            }
            u_ = u;
        }

        Eigen::VectorXd& getPredict(void) { return x_pred_; };
        Eigen::VectorXd& getEstimate(void) { return x_est_; };

        /** @brief Performs the KF prediction step
         *
         * Calculates the state predicted from the model and updates the estimate
         * error covariance matrix.
         */
        void predict(void) {
            x_pred_ = A_ * x_est_ + B_ * u_;
            P_ = A_ * P_ * A_.transpose() + Q_;
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
            K_ = P_ * H_.transpose() * (H_ * P_ * H_.transpose() + R_).inverse();
            x_est_ = x_pred_ + K_ * (z - H_ * x_pred_);
            P_ = (I_ - K_ * H_) * P_;
        }

    private:
        // Matrices for computation
        Eigen::MatrixXd A_, B_, H_, Q_, R_;

        size_t l_; // dimension of input vector. Can be zero for autonomous system.
        size_t m_; // dimension of output vector. Can be zero for prediction only.
        size_t n_; // dimension of state vector

        // n-size identity
        Eigen::MatrixXd I_;
        Eigen::MatrixXd P_, K_;

        // input vector
        Eigen::VectorXd u_;

        // Estimated states
        Eigen::VectorXd x_pred_, x_est_;
};

}

#endif  //!__KALMAN_FILTER__H__