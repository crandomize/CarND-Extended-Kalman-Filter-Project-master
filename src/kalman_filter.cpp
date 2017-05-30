#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si; 
	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // To calculate z_pred we need to convert from cartesian to polar.
  double rho_pred, phi_pred, rhodot_pred;
  rho_pred = sqrt(pow(x_[0], 2) + pow(x_[1], 2));

  phi_pred = 0.0;
  if (fabs(x_[0]) > 0.001) {
    phi_pred  = atan2(x_[1], x_[0]);
  }

  rhodot_pred = 0.0;
  if (fabs(rho_pred) > 0.001) {
    rhodot_pred = (x_[0] * x_[2] + x_[1] * x_[3]) / rho_pred;
  }
  VectorXd z_pred(3);
  z_pred << rho_pred, phi_pred, rhodot_pred;

	//VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;

  //if radar measurement (size is 3), normalize angle
  if (y.size() == 3) y(1) = atan2(sin(y(1)), cos(y(1))); 

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
