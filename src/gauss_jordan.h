#ifndef GAUSS_JORDAN_H
#define GAUSS_JORDAN_H

#include <Eigen/Dense>

bool check_if_inversible (const Eigen::MatrixXf&);

void compute_inverse (Eigen::MatrixXf&);

#endif
