#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <deque>
#include <set>
#include <string>
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

void print_peak(const Eigen::MatrixXf& A, const std::vector<int>& path)
{
    std::cout << "peak = " << peak_mem(A, path) << std::endl;
}

void print_peak_comparison(const Eigen::MatrixXf& A, const std::vector<int>& path1,
                           const std::vector<int>& path2)
{
    std::cout << "peak = " << peak_mem(A, path1) << std::endl;
    std::cout << "peak (random path) = " << peak_mem(A, path2) << std::endl;
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

/*******************************************************************************
 * Lists each integer of a particular string.
 *
 *
 * @param s String.
 * @param del Delimiter between integers.
 ******************************************************************************/
std::vector<int> string_split(const std::string& s, const char del)
{
    std::vector<int> tokens;
    std::string token;
    std::istringstream token_stream(s);

    while(std::getline(token_stream, token, del))
        {
            tokens.push_back(std::stoi(token));
        }

    return tokens;
}

/*******************************************************************************
 * From the input file, generates a matrix A representing the graph and a vector
 * O in order to check the ordering when applying the chosen algorithm.
 *
 * Input file has the following format:
 * First line corresponds to the number of nodes (N)
 * From the second line onwards:
 * n1 o1 (i.e., node 1 has output size o1)
 * n2 o2 n1 (i.e., node 2 has output size o2 and receives n1 as input)
 *
 * @param file_name String containing file's name.
 * @param A Adjacency matrix of the graph.
 * @param O Ordering constraints (given in the form L_x > [L_yi]).
 ******************************************************************************/
void read_file(const std::string& file_name, Eigen::MatrixXf& A, int N)
{
    std::vector<std::vector<int>> tokens(N), O(N);

    // File initialization
    std::ifstream file(file_name);

    if (!file.is_open()) return;

    // Read first line and ignore (we already know N)
    std::string line;
    std::getline(file, line);

    // Reading file line by line
    int i = 0;

    while (std::getline(file, line))
    {
        tokens[i] = string_split(line, ' ');

        A(i, i) = float(tokens[i][0]);

        for (std::vector<int>::size_type j = 1; j < tokens[i].size(); ++j)
        {
            O[i].push_back(tokens[i][j]);
        }

        i++;
    }

    // Update matrix A only after all data is read
    for (std::vector<std::vector<int>>::size_type i = 0; i < tokens.size(); ++i)
    {
        for (std::vector<int>::size_type j = 1; j < tokens[i].size(); ++j)
        {
            A(i, tokens[i][j]) = A(tokens[i][j], tokens[i][j]);
        }
    }

    file.close();
}

int read_N_from_file (const std::string& file_name)
{
    std::ifstream file(file_name);
    if (!file.is_open()) return -1;

    std::string line;
    std::getline(file, line);
    int N = std::stoi(line);

    file.close();

    return N;
}
