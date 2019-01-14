//Lara Querciagrossa
#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <vector>

double* crossprod(double crd1[], double crd2[])
{
  //Compute cross product between two vectors (double[3])

  double* cross = new double[3];
  for (int i=0; i<3; i++)
    {
      cross[i]=0.0;
    }

  cross[0] = crd1[1]*crd2[2]-crd1[2]*crd2[1];
  cross[1] = crd1[2]*crd2[0]-crd1[0]*crd2[2];
  cross[2] = crd1[0]*crd2[1]-crd1[1]*crd2[0];

  return cross;
}


double* MatrixToQuaternion (double matrix[3][3], bool active)
{
  std::vector<double> quat;

  //Convert a rotation matrix to the corresponding quaternion

  double n0 = sqrt(matrix[0][0]*matrix[0][0] + matrix[0][1]*matrix[0][1] + matrix[0][2]*matrix[0][2]);
  double n1 = sqrt(matrix[1][0]*matrix[1][0] + matrix[1][1]*matrix[1][1] + matrix[1][2]*matrix[1][2]);
  double n2 = sqrt(matrix[2][0]*matrix[2][0] + matrix[2][1]*matrix[2][1] + matrix[2][2]*matrix[2][2]);
  matrix[0][0] = matrix[0][0]/n0;
  matrix[0][1] = matrix[0][1]/n0;
  matrix[0][2] = matrix[0][2]/n0;
  matrix[1][0] = matrix[1][0]/n1;
  matrix[1][1] = matrix[1][1]/n1;
  matrix[1][2] = matrix[1][2]/n1;
  matrix[2][0] = matrix[2][0]/n2;
  matrix[2][1] = matrix[2][1]/n2;
  matrix[2][2] = matrix[2][2]/n2;

  double trace = matrix[0][0]+matrix[1][1]+matrix[2][2];
  double qw, qx, qy, qz;

  if (trace > 0.0)
    {
      qw = sqrt(1+trace)/2.0;
      qx = (matrix[2][1] - matrix[1][2])/(4*qw);
      qy = (matrix[0][2] - matrix[2][0])/(4*qw);
      qz = (matrix[1][0] - matrix[0][1])/(4*qw);
    }
  else if (matrix[0][0] > matrix[1][1] && matrix[0][0] > matrix[2][2])
    {
      double tmp = matrix[0][0] - matrix[1][1] - matrix[2][2];
      qx = sqrt(1.0+tmp)/2.0;
      qw = (matrix[2][1] - matrix[1][2])/(4*qx);
      qy = (matrix[0][1] + matrix[1][0])/(4*qx);
      qz = (matrix[0][2] + matrix[2][0])/(4*qx);
    }
  else if (matrix[1][1] > matrix[2][2])
    {
      double tmp = matrix[1][1] - matrix[0][0] - matrix[2][2];
      qy = sqrt(1.0+tmp)/2.0;
      qw = (matrix[0][2] - matrix[2][0])/(4*qy);
      qx = (matrix[0][1] + matrix[1][0])/(4*qy);
      qz = (matrix[1][2] + matrix[2][1])/(4*qy);
    }
  else 
    {
      double tmp = matrix[2][2] - matrix[0][0] - matrix[1][1];
      qz = sqrt(1.0+tmp)/2.0;
      qw = (matrix[1][0] - matrix[0][1])/(4*qz);
      qx = (matrix[0][2] + matrix[2][0])/(4*qz);
      qy = (matrix[1][2] + matrix[2][1])/(4*qz);
    }

  if (active == 0)
    {
      qw = -qw;
      qx = qx;
      qy = qy;
      qz = qz;
    }

  quat.push_back(qw);
  quat.push_back(qx);
  quat.push_back(qy);
  quat.push_back(qz);

  double* q = new double[4];
  for (int i=0; i<4; i++)
    {
      q[i] = 0.0;
    }

  q[0] = quat[0];
  q[1] = quat[1];
  q[2] = quat[2];
  q[3] = quat[3];

  return q;
}


