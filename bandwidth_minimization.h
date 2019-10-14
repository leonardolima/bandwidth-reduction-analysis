#ifndef GAUSS_JORDAN_H
#define GAUSS_JORDAN_H

#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <chrono>
#include <Eigen/Dense>
#include "gauss_jordan.h"

bool node_has_label (const Eigen::MatrixXf&, unsigned int);

std::vector<std::pair<int, int>> sort_row_deg (const Eigen::VectorXf&, 
                                               const std::vector<int>&, 
                                               const Eigen::MatrixXf&);

int node_index (const Eigen::MatrixXf&, int);

void nodal_numbering (const Eigen::MatrixXf&, Eigen::MatrixXf&, const std::vector<int>&);

int matrix_deg (Eigen::MatrixXf&);

void compareMatrices(const Eigen::MatrixXf&, const Eigen::MatrixXf&);

std::vector<int> compute_row_deg (const Eigen::MatrixXf&);

std::vector<int> select_starting_nodes (Eigen::MatrixXf&);

void compute_results (std::vector<int>&, const Eigen::MatrixXf&,
                      Eigen::MatrixXf&, Eigen::MatrixXf&, Eigen::MatrixXf&);

void generate_binary_random_matrix (Eigen::MatrixXf&, int);

#endif
