#include "kalman_filter.h"
#define PI 3.14159265

#include <iostream>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init( VectorXd &x_in, 
    MatrixXd &P_in, 
    MatrixXd &F_in,
    MatrixXd &H_in, 
    MatrixXd &Hj_in, 
    MatrixXd &R_in, 
    MatrixXd &R_ekf_in, 
    MatrixXd &Q_in ) 
{
  cout << "In KalmanFilter::Init" << endl;
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;

  R_ = R_in;



  Q_ = Q_in;

  //define an identity matrix


}

void KalmanFilter::Predict() 
{
  // predict the state
  x_ = F_*x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_*P_*Ft + Q_;
}


void KalmanFilter::Update(const VectorXd &z) 
{
  // Update the state using Kalman Filter equations

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Ht*Si;

  // New state
  x_ = x_ + ( K*y );
  P_ -= K * H_ * P_;

}



void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];


//avoid division by zero
  if( px == 0. && py == 0. )
    return;

  // Hj_ = tools.CalculateJacobian( x_ );


  VectorXd z_pred(3);

  double rho = sqrt( px*px + py*py );
  double phi = atan2( py, px ); 
  double rho_dot;


  rho_dot = (px *vx + py*vy)/rho;



  z_pred << rho, phi, rho_dot;


  // Update the state using Extended Kalman Filter equations
  VectorXd y = z - z_pred;
  if( y[1] > PI )
    y[1] -= 2.f*PI;
  if( y[1] < -PI )
    y[1] += 2.f*PI;


  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Ht*Si;

  // Compute new state
  x_ = x_ + ( K*y );
  P_ -= K * H_ * P_;
}



