#ifndef BANDWIDTH_MINIMIZATION_H
#define BANDWIDTH_MINIMIZATION_H

#include <vector>
#include <Eigen/Dense>

bool node_has_label (const Eigen::MatrixXf&, unsigned int);

std::vector<std::pair<int, int>> sort_row_deg (const Eigen::VectorXf&, 
                                               const std::vector<int>&, 
                                               const Eigen::MatrixXf&);

int node_index (const Eigen::MatrixXf&, int);

void nodal_numbering (const Eigen::MatrixXf&, Eigen::MatrixXf&, const std::vector<int>&);

int matrix_deg (Eigen::MatrixXf&);

void compare_matrices(const Eigen::MatrixXf&, const Eigen::MatrixXf&);

std::vector<int> compute_rows_deg (const Eigen::MatrixXf&);

std::vector<int> select_starting_nodes (Eigen::MatrixXf&);

void compute_matrices (std::vector<int>&, Eigen::MatrixXf&,
                      Eigen::MatrixXf&, Eigen::MatrixXf&);

void generate_binary_random_matrix (Eigen::MatrixXf&);

void run_tests (int);

void run_algorithm (const Eigen::MatrixXf&, Eigen::MatrixXf&, Eigen::MatrixXf&);

#endif
