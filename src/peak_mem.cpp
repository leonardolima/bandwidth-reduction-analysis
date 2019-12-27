/******************************
 *                            *
 *    Leonardo Lima, 2019     *
 *                            *
/******************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "cuthill_mckee.h"

/* Function: string_split
 *
 * Inputs:
 *         s   - string
 *         del - delimiter (in our particular case it is a single space)
 *
 * List each integer on the string s
 *
 */

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

/* Function: read_file
 *
 * Inputs:
 *
 * Reads a .txt file in the following format:
 * First line corresponds to the number of nodes (N)
 * From the second line onwards:
 * n1 - o1 (i.e., node 1 has output size o1)
 * n2 - o2 n1 (i.e., node 2 has output size o2 and receives n1 as input)
 *
 * From there, generates a matrix A representing the graph and a vector O
 * in order to check the ordering when applying the algorithm
 *
 */
void read_file (const std::string& file_name, Eigen::MatrixXf& A, std::vector<std::vector<int>>& O)
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
        std::vector tokens = string_split(line, ' ');

        A(i, i) = tokens[0];

        for (std::vector<int>::size_type j = 1; j < tokens.size(); ++j)
        {
            O[i].push_back(tokens[j]);
            A(tokens[j], i) = A(tokens[j], tokens[j]);
            A(i, tokens[j]) = A(tokens[j], tokens[j]);
        }

        i++;
    }

    file.close();

    std::cout << "A = " << std::endl;
    std::cout << A << std::endl;
}

void peak_mem(const std::string& file_name)
{
    int N = read_N_from_file(file_name);

    // Adjacency matrix representing the graph (symmetric)
    Eigen::MatrixXf A = Eigen::MatrixXf::Identity(N, N);

    // Vector corresponding to the ordering
    std::vector<std::vector<int>> O(N);

    read_file(file_name, A, O);
}
