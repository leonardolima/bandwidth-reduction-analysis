#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <Eigen/Dense>

void print_set (const std::set<int>& s)
{
    for (auto it = s.begin(); it != s.end(); ++it)
    {
        std::cout << ' ' << *it;
    }
    std::cout << std::endl;
}

void print_queue (const std::queue<int>& q)
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

void print_map (const std::map<std::pair<int, int>, int>& m)
{
    for (const auto &it : m)
    {
        auto& k = it.first;
        auto& v = it.second;
        std::cout << "m[" << k.first << ", " << k.second  << "] = (" << v << ") " << std::endl;
    }
}

void matrix_to_csv (const Eigen::MatrixXf& R)
{
    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    // Initialize output file
    std::ofstream f("out.csv");

    if (!f.is_open()) return;

    f << R.format(CSVFormat);

    f.close();
}
