#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  ///* time when the state is true, in us
  time_us_ = 0;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  ///* Predicted sigma points dimension
  n_sig_ = 2 * n_aug_ + 1;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  //set vector for weights
  weights_ = VectorXd(n_sig_);
  weights_.setConstant(n_sig_, (0.5/(lambda_+n_aug_)));
  weights_(0) = lambda_ / (lambda_+n_aug_);
  ///* set measurement dimension, lidar can measure px and py
  n_z_laser_ = 2;

  ///* define covariance matrix R_laser_ for lidar measurement
  R_laser_ = MatrixXd(n_z_laser_,n_z_laser_);
  R_laser_.fill(0.0);
  R_laser_(0,0) = std_laspx_ * std_laspx_;
  R_laser_(1,1) = std_laspy_ * std_laspy_;

  ///* set measurement dimension, radar can measure r, phi, and r_dot
  n_z_radar_ = 3;

  ///* define covariance matrix R_radar_ for radar measurement
  R_radar_ = MatrixXd(n_z_radar_,n_z_radar_);
  R_radar_.fill(0.0);
  R_radar_(0,0) = std_radr_ * std_radr_;
  R_radar_(1,1) = std_radphi_ * std_radphi_;
  R_radar_(2,2) = std_radrd_ * std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  // define constants
  double delta_t;

  if(!is_initialized_){
    /**
    TODO:
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    // first measurement
    cout<<"UKF: "<<endl;
    // store Normalized Innovation Squared (NIS)
    ofstream output_radar, output_laser;
    output_radar.open("NIS_RADAR.txt");
    output_radar << "NIS_RADAR" << "\n";
    output_radar.close();

    output_laser.open("NIS_LASER.txt");
    output_laser << "NIS_LASER" << "\n";
    output_laser.close();

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rhod = meas_package.raw_measurements_[2];

      x_ << rho * cos(phi), rho * sin(phi), rhod, 0, 0;

      //cout<<"RADAR: "<<x_<<endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      //cout<<"LASER: "<<x_<<endl;
    }

    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
  }
  else{
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_){
      delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
      time_us_ = meas_package.timestamp_;
      // Prediction
      Prediction(delta_t);

      // Update
      UpdateRadar(meas_package);
    }
    else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_){
      delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
      time_us_ = meas_package.timestamp_;
      // Prediction
      Prediction(delta_t);

      // Update
      UpdateLidar(meas_package);
    }
  }

  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // define temp constants
  double tmp_sqr = sqrt(lambda_+n_aug_);
  double half_delta_t_sq = 0.5 * delta_t * delta_t;
  VectorXd xsig(n_aug_), xsig1(n_x_), xsig2(n_x_), x_diff(n_x_);
  double tmp1, x2ox4;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  //calculate square root of P_
  MatrixXd A = P_.llt().matrixL();

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd SRM = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i=0;i<n_aug_;i++){
      Xsig_aug.col(i+1) = x_aug + tmp_sqr * SRM.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - tmp_sqr * SRM.col(i);
  }

  for(int i=0;i<n_sig_;i++){
    //predict sigma points
    xsig = Xsig_aug.col(i);

    //avoid division by zero
    if(fabs(xsig(4))<0.001){
        xsig1 <<
        xsig(2)*cos(xsig(3))*delta_t,
        xsig(2)*sin(xsig(3))*delta_t,
        0,
        0,
        0;

        Xsig_pred_.col(i) = xsig.head(5) + xsig1;
    }
    else{
        tmp1 = xsig(3)+xsig(4)*delta_t;
        x2ox4 = xsig(2)/xsig(4);

        xsig1 <<
        x2ox4*(sin(tmp1) - sin(xsig(3))),
        x2ox4*(-cos(tmp1) + cos(xsig(3))),
        0,
        xsig(4)*delta_t,
        0;

        Xsig_pred_.col(i) = xsig.head(5) + xsig1;
    }
    xsig2 <<
        half_delta_t_sq*cos(xsig(3))*xsig(5),
        half_delta_t_sq*sin(xsig(3))*xsig(5),
        delta_t*xsig(5),
        half_delta_t_sq*xsig(6),
        delta_t*xsig(6);
    Xsig_pred_.col(i) += xsig2;
  }

  //predict state mean
  x_ = Xsig_pred_ * weights_;

  //predict state covariance matrix
  P_.fill(0.0);
  for(int i=0;i<n_sig_;i++){
      x_diff = Xsig_pred_.col(i) - x_;

      //normalize angle
      NormalizeAngle(x_diff(3));

      P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Laser updates
  cout<<"\nLASER updates\n";

  // define temp constants
  VectorXd z_diff(n_z_laser_), z(n_z_laser_), x_diff(n_x_);
  double px, py;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z_laser_, n_sig_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_laser_);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_laser_, n_z_laser_);

  //calculate mean predicted measurement
  z_pred = Zsig * weights_;

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for(int i=0;i<n_sig_;i++){
      z_diff = Zsig.col(i) - z_pred;

      S += weights_(i) * z_diff * z_diff.transpose();
  }

  S += R_laser_;

  //create example vector for incoming radar measurement
  z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0;i < n_sig_;i++){
      x_diff = Xsig_pred_.col(i) - x_;
      z_diff = Zsig.col(i) - z_pred;

      Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  z_diff = z - z_pred;

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

  // calculate Normalized Innovation Squared (NIS)
  float NIS_laser = z_diff.transpose() * S.inverse() * z_diff;

  // store Normalized Innovation Squared (NIS)
  ofstream output_laser;
  output_laser.open("NIS_LASER.txt", ios_base::app);
  // output the NIS values
  output_laser << NIS_laser << "\n";
  output_laser.close();

  return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Radar updates
  cout<<"\nRADAR updates\n";

  // define temp constants
  VectorXd z_diff(n_z_radar_), z(n_z_radar_), x_diff(n_x_);
  double px, py, nu, psi, rho, phi, rhod;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, n_sig_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);

  //transform sigma points into measurement space
  for(int i=0;i<n_sig_;i++){
      px = Xsig_pred_(0,i);
      py = Xsig_pred_(1,i);
      nu = Xsig_pred_(2,i);
      psi = Xsig_pred_(3,i);

      rho = sqrt(px*px + py*py);

      // define a special case for atan2(0, 0)
      if((px == 0) && (py == 0)){
        phi = 0;
      }
      else{
        phi = atan2(py, px);
      }

      // avoid dividing by 0
      if(rho < 0.0001){
        rhod = (px*cos(psi)* nu + py*sin(psi)* nu)  / 0.0001;
      }
      else{
        rhod = (px*cos(psi)* nu + py*sin(psi)* nu)  / rho;
      }

      Zsig.col(i) << rho, phi, rhod;
  }

  //calculate mean predicted measurement
  z_pred = Zsig * weights_;

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for(int i=0;i<n_sig_;i++){
      z_diff = Zsig.col(i) - z_pred;

      NormalizeAngle(z_diff(2));

      S += weights_(i) * z_diff * z_diff.transpose();
  }

  S += R_radar_;

  //create example vector for incoming radar measurement
  z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0;i < n_sig_;i++){
      x_diff = Xsig_pred_.col(i) - x_;
      z_diff = Zsig.col(i) - z_pred;

      NormalizeAngle(x_diff(3));
      NormalizeAngle(z_diff(1));

      Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  z_diff = z - z_pred;

  NormalizeAngle(z_diff(1));

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // print the output
  cout << "x_ = \n" << x_ << endl;
  cout << "P_ = \n" << P_ << endl;

  // calculate Normalized Innovation Squared (NIS)
  float NIS_laser = z_diff.transpose() * S.inverse() * z_diff;

  // store Normalized Innovation Squared (NIS)
  ofstream output_radar;
  output_radar.open("NIS_RADAR.txt", ios_base::app);
  // output the NIS values
  output_radar << NIS_laser << "\n";
  output_radar.close();

  return;
}

/**
 * Normalize angle.
 */
void UKF::NormalizeAngle(double &angle){
  if((angle > M_PI) || (angle <= -M_PI)){
    angle -= (2*M_PI) * floor((angle + M_PI) / (2*M_PI));
  }
}
