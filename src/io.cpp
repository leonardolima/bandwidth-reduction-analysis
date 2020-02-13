#include <fstream>
#include <iostream>
#include <map>
#include <Eigen/Dense>

void print_map(const std::map<std::pair<int, int>, int>& m)
{
    for (const auto &it : m)
    {
        auto& k = it.first;
        auto& v = it.second;
        std::cout << "m[" << k.first << ", " << k.second  << "] = (" << v << ") " << std::endl;
    }
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
