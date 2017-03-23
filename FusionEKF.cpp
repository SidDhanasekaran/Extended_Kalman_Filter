#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// Constructor
FusionEKF::FusionEKF() {

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // Initialize matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

}

// Destructor
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // Initialize FusionEKF state variables
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Initial transition matrix, F
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    // State covariance matrix, P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    // Initial Q matrix
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               1, 0, 1, 0,
               0, 1, 0, 1;

    noise_ax = 9;
    noise_ay = 9;

    // Laser measurement matrix
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    // Radar measurement matrix
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

    Hj_ << 1, 1, 0, 0,
           1, 1, 0, 0,
           1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      // Rho = Range
      float rho = measurement_pack.raw_measurements_[0];

      // Phi = Bearing
      float phi = measurement_pack.raw_measurements_[1];

      //Adding range rate

      float range_rate = measurement_pack.raw_measurements_[2];

      float cartesian_x = rho * cos(phi);
      float cartesian_y = rho * sin(phi);
      float vx = range_rate *cos(phi);
      float vy = range_rate *sin(phi);

      if (cartesian_x == 0 or cartesian_y == 0){
        return;
      }

      ekf_.x_ << cartesian_x, cartesian_y, vx, vy;

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

        if (measurement_pack.raw_measurements_[0] == 0 or measurement_pack.raw_measurements_[1] == 0){
          return;
        }

        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

    }
    is_initialized_ = true;

    return;
    }

    // Single measurement

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    if (dt > 0.001) {
      // Update state transistion matrix with new time
      ekf_.F_ << 1, 0, dt, 0,
                 0, 1, 0, dt,
                 0, 0, 1, 0,
                 0, 0, 0, 1;

      ekf_.Q_ << (pow(dt, 4) / 4 * noise_ax), 0, (pow(dt, 3) / 2 * noise_ax), 0,
                  0, (pow(dt, 4) / 4 * noise_ay), 0, (pow(dt, 3) / 2 * noise_ay),
                  (pow(dt, 3) / 2 * noise_ax), 0, pow(dt, 2) * noise_ax, 0,
                  0, (pow(dt, 3) / 2 * noise_ay), 0, pow(dt, 2) * noise_ay;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

      ekf_.Predict();
    }
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Use the sensor type to perform the update step.
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;

    // Rho = Range
    double rho = sqrt(pow(ekf_.x_[0], 2) + pow(ekf_.x_[1], 2));

    // Phi = Bearing
    double phi = atan(ekf_.x_[1] / ekf_.x_[0]);

    // Rho_dot = Range rate
    double rho_dot = ((ekf_.x_[0] * ekf_.x_[2] + ekf_.x_[1] * ekf_.x_[3]) / (sqrt(pow(ekf_.x_[0], 2) + pow(ekf_.x_[1], 2))));

    MatrixXd zpred(3, 1);
    zpred << rho, phi, rho_dot;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_, zpred);

  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;

    // Update Laser state and covariance matrices
    ekf_.Update(measurement_pack.raw_measurements_);

  }

}
