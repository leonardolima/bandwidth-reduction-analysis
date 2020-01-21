#ifndef CUTHILL_MCKEE_H
#define CUTHILL_MCKEE_H

#include <Eigen/Dense>

void compute_matrices (const Eigen::MatrixXf&, Eigen::MatrixXf&, Eigen::MatrixXf&, const std::vector<std::vector<int>>&);

#endif
