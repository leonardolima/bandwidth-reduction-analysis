#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <deque>
#include <set>
#include <Eigen/Dense>
#include "io.h"
#include "basic.h"
#include "peak_mem.h"

void print_set(const std::set<int>& s)
{
    for (auto it = s.begin(); it != s.end(); ++it)
    {
        std::cout << ' ' << *it;
    }
    std::cout << std::endl;
}

void print_queue(const std::queue<int>& q)
{
    std::queue<int> q_copy = q;

    std::cout << "q = ";
    while (!q_copy.empty())
    {
        std::cout << q_copy.front() << " ";
        q_copy.pop();
    }
    std::cout << std::endl;
}

void print_deque(const std::deque<int>& d)
{
    for (std::deque<int>::size_type i = 0; i < d.size(); ++i)
    {
        std::cout << "d[" << i << "] = " << d[i] << std::endl;
    }
}

void print_map(const std::map<std::pair<int, int>, int>& m)
{
    for (const auto &it : m)
    {
        auto& k = it.first;
        auto& v = it.second;
        std::cout << "m[" << k.first << ", " << k.second  << "] = (" << v << ") " << std::endl;
    }
}

void print_bandwidth_comparison(const Eigen::MatrixXf& A, const Eigen::MatrixXf& R)
{
    std::cout << "dim(A) = dim(R) = " << A.rows() << "x" << A.cols() << std::endl;
    std::cout << "bandwidth(A) = " << matrix_bandwidth(A) << std::endl;
    std::cout << "bandwidth(R) = " << matrix_bandwidth(R) << std::endl;
}

void print_peak_comparison(const Eigen::MatrixXf& A, const std::vector<int>& path1,
                           const std::vector<int>& path2)
{
    std::cout << "peak_1 = " << peak_mem(A, path1) << std::endl;
    std::cout << "peak_2 = " << peak_mem(A, path2) << std::endl;
}

void matrix_to_csv(const Eigen::MatrixXf& R)
{
    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    // Initialize output file
    std::ofstream f("out.csv");

    if (!f.is_open()) return;

    f << R.format(CSVFormat);

    f.close();
}
