#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include <Eigen/Dense>

void ludcmp(Eigen::MatrixXf&, Eigen::VectorXf&, double&);

void lubksb(const Eigen::MatrixXf&, const Eigen::VectorXf&, Eigen::VectorXf&);

#endif
