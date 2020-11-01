#ifndef BASIC_H
#define BASIC_H

#include <Eigen/Dense>

int matrix_bandwidth (const Eigen::MatrixXf&);

void parents_from_matrix(const Eigen::MatrixXf&, std::vector<std::vector<int>>&);

void children_from_matrix(const Eigen::MatrixXf&, std::vector<std::vector<int>>&);

void mem_from_matrix(const Eigen::MatrixXf&, std::vector<int>&);

void apply_symmetry (Eigen::MatrixXf&);

#endif
