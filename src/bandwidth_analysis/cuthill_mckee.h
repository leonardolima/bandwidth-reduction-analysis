#ifndef BANDWIDTH_MINIMIZATION_H
#define BANDWIDTH_MINIMIZATION_H

#include <vector>
#include <Eigen/Dense>

void run_tests (int);

void run_algorithm (const Eigen::MatrixXf&, Eigen::MatrixXf&, Eigen::MatrixXf&);

bool node_has_label (const Eigen::MatrixXf&, unsigned int);

int node_label (const Eigen::MatrixXf&, unsigned int);

int node_index (const Eigen::MatrixXf&, int);

std::vector<std::pair<int, int>> sort_rows_deg (const Eigen::VectorXf&,
                                                const std::vector<int>&,
                                                const Eigen::MatrixXf&);

int matrix_deg (const Eigen::MatrixXf&);

std::vector<int> compute_rows_deg (const Eigen::MatrixXf&);

#endif
