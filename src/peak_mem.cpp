#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "adapted_cuthill_mckee.h"

/*******************************************************************************
 * Lists each integer of a particular string.
 *
 *
 * @param s String.
 * @param del Delimiter between integers.
 ******************************************************************************/
std::vector<int> string_split (const std::string& s, const char del)
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
 * @param A Adjacency matrix of the graph (symmetric).
 * @param O Ordering constraints (given in the form L_x > [L_yi]).
 ******************************************************************************/
void read_file (const std::string& file_name, Eigen::MatrixXf& A,
                std::vector<std::vector<int>>& O)
{
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
        std::vector<int> tokens = string_split(line, ' ');

        A(i, i) = float(tokens[0]);

        for (std::vector<int>::size_type j = 1; j < tokens.size(); ++j)
        {
            O[i].push_back(tokens[j]);
            A(tokens[j], i) = A(tokens[j], tokens[j]);
            A(i, tokens[j]) = A(tokens[j], tokens[j]);
        }

        i++;
    }

    file.close();
}

/*******************************************************************************
 * Changes constraints format.
 *
 * Constraints are given in the form L_x > [L_yi], but in order to keep search
 * O(1) the format is changed to L_x < [L_yi].
 *
 *
 * @param O     Ordering constraints (in the form L_x > [L_yi]).
 * @param new_O Ordering constraints (in the form L_x < [L_yi]).
 ******************************************************************************/
void convert_vector(const std::vector<std::vector<int>>& O,
                    std::vector<std::vector<int>>& new_O)
{
    for (std::vector<std::vector<int>>::size_type i = 0; i < O.size(); ++i)
    {
        for (std::vector<int>::size_type j = 0; j < O[i].size(); ++j)
        {
            new_O[O[i][j]].push_back(i);
        }
    }
}

void peak_mem(const std::string& file_name)
{
    int N = read_N_from_file(file_name);

    // Adjacency matrix representing the graph (symmetric)
    Eigen::MatrixXf A = Eigen::MatrixXf::Identity(N, N);

    // Permutation matrix (in case of Cuthill-McKee algorithm)
    Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N, N);

    // Resulting matrix
    Eigen::MatrixXf R = Eigen::MatrixXf::Zero(N, N);

    // Vector corresponding to the ordering L_x > [L_yi]
    std::vector<std::vector<int>> O(N);

    // Vector corresponding to the ordering L_x < [L_yi]
    std::vector<std::vector<int>> new_O(N);

    read_file(file_name, A, O);

    convert_vector(O, new_O);

    apply_adapted_cuthill_mckee(A, P, R, new_O);
}
