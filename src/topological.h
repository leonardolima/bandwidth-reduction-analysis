#ifndef TOPOLOGICAL_H
#define TOPOLOGICAL_H

#include <Eigen/Dense>

void apply_topological(const Eigen::MatrixXf&, Eigen::MatrixXf&,
                       std::vector<int>&, std::vector<int>&);

#endif
