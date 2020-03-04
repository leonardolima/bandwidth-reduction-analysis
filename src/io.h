#ifndef IO_H
#define IO_H

#include <Eigen/Dense>

void matrix_to_csv(const Eigen::MatrixXf&);

void print_queue(const std::queue<int>&);

void print_set (const std::set<int>& s);

#endif
