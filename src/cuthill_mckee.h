#ifndef CUTHILL_MCKEE_H
#define CUTHILL_MCKEE_H

#include <Eigen/Dense>

void run_algorithm (const Eigen::MatrixXf&, Eigen::MatrixXf&, Eigen::MatrixXf&, std::vector<std::vector<int>>&);

#endif
