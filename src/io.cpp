/******************************
 *                            *
 *    Leonardo Lima, 2019     *  
 *                            *
/******************************/

#include <fstream>
#include <Eigen/Dense>

/* Function: matrix_to_csv
 *
 * Inputs: R - matrix
 *
 * The matrix is saved as a CSV file
 */
void matrix_to_csv(const Eigen::MatrixXf& R)
{
    const Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    // Initialize output file
    std::ofstream f("out.csv");

    if (!f.is_open()) return;

    f << R.format(CSVFormat);

    f.close();
}
