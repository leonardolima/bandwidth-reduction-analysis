#ifndef LEVEL_H
#define LEVEL_H

#include <Eigen/Dense>

void apply_levels(const Eigen::MatrixXf&, Eigen::MatrixXf&,
                  std::vector<int>&);

#endif
