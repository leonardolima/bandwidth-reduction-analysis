#ifndef CUTHILL_MCKEE_H
#define CUTHILL_MCKEE_H

#include <vector>
#include <Eigen/Dense>

void run_tests(int);

void run_algorithm(const Eigen::MatrixXf&, Eigen::MatrixXf&, Eigen::MatrixXf&);

bool node_has_label(const Eigen::MatrixXf&, unsigned int);

int node_label(const Eigen::MatrixXf&, unsigned int);

int node_index_adj_matrix(const Eigen::MatrixXf&, int);

std::vector<std::pair<int, int>> sort_adj_nodes(const Eigen::VectorXf&,
                                                const std::vector<int>&,
                                                const Eigen::MatrixXf&);

int matrix_deg(const Eigen::MatrixXf&);

#endif
