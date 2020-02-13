#ifndef CUTHILL_MCKEE_H
#define CUTHILL_MCKEE_H

#include <Eigen/Dense>

void apply_adapted_cuthill_mckee (const Eigen::MatrixXf&, Eigen::MatrixXf&, Eigen::MatrixXf&, const std::vector<std::vector<int>>&);

#endif
