#ifndef TOPOLOGICAL_H
#define TOPOLOGICAL_H

#include <Eigen/Dense>

void apply_topological_all(const Eigen::MatrixXf&, Eigen::MatrixXf&,
                           std::vector<int>&);

void apply_topological(const Eigen::MatrixXf&, std::vector<int>&);

#endif
