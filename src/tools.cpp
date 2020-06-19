#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() : mse(4), 
    lastmse(4),
    residual(4),
    kahanerror(4),
    rmse(4)
{
  resetRMSE();
}

Tools::~Tools() {}

void Tools::resetRMSE()
{
  mse << 0, 0, 0, 0;
  lastmse << 0, 0, 0, 0;
  residual << 0, 0, 0, 0;
  kahanerror << 0, 0, 0, 0;
}

// Calculate the RMSE
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  // Calculate a Jacobian for the transformation from the state vector 
  // px, py, vx, vy to the radar measurement space
  // rho, phi, rhodot.
  

  // This should be copy-elided into wherever it's being returned.
  MatrixXd Hj(3,4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  if( px == 0 && py == 0 )
  {
    cout << "Error:  division by zero in CalculateJacobian" << endl;
    return Hj;
  }

  //compute the Jacobian 
  float rho = sqrt( px*px + py* py );
  float rho2 = rho*rho;
  float rho3 = rho2*rho;
  Hj <<                 px/rho,                    py/rho,      0,      0,
                      -py/rho2,                   px/rho2,      0,      0,
     py*( vx*py - vy*px )/rho3, px*( vy*px - vx*py )/rho3, px/rho, py/rho;

  return Hj;
}
