#ifndef BASIC_H
#define BASIC_H

#include <Eigen/Dense>

int matrix_bandwidth (const Eigen::MatrixXf&);

void parents_from_matrix(const Eigen::MatrixXf&, std::vector<std::vector<int>>&);

void children_from_matrix(const Eigen::MatrixXf&, std::vector<std::vector<int>>&);

void mem_from_matrix(const Eigen::MatrixXf&, std::vector<int>&);

void apply_symmetry (Eigen::MatrixXf&);

std::vector<int> compute_nodes_deg(const Eigen::MatrixXf&);

void label_sorted_nodes(const Eigen::MatrixXf&, Eigen::MatrixXf&,
                        std::vector<std::pair<int, int>>);

std::vector<std::pair<int, int>> make_sorted_pairs(const std::vector<int>);

#endif
