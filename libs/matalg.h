#ifndef MATALG_H
#define MATALG_H

#include <vector>
#include <stdexcept>
#include <cmath>

/*
Matrices are defined as std::vector<std::vector<double>>
Second containers are the rows. i.e., In std::vector<std::vector<double>> A, for A[i][j],
'i' is the row number, 'j' is the column number.
*/

//Validity functions
namespace matalg
{

bool validate_matrix(const std::vector<std::vector<double>>&);
//Returns true if every row has same number of elements
bool validate_squareMat(const std::vector<std::vector<double>>&);
//Returns true if the matrix is square.
//If this is called, calling validate_matrix is not required.
bool validate_equal_dimension(const std::vector<std::vector<double>>&,
                               const std::vector<std::vector<double>>&);
//Returns true if two matrices have same dimensions.
//Should be used with validate_matrix.

//Operation functions

//Linear system
//Determinant
double det(const std::vector<std::vector<double>>&, bool validate = true);
//Linear system solution
std::vector<double> solve_linSystem(const std::vector<std::vector<double>>&, const std::vector<double>&);
//Solves a linear system in the form Ax=B

//Algebraic
//A+B
std::vector<std::vector<double>> mat_add(const std::vector<std::vector<double>>&,
                                             const std::vector<std::vector<double>>&);
//A-B
std::vector<std::vector<double>> mat_sub(const std::vector<std::vector<double>>&,
                                            const std::vector<std::vector<double>>&);
//mA
std::vector<std::vector<double>> mat_scale(const std::vector<std::vector<double>>&, double);
//A*B
std::vector<std::vector<double>> mat_mul(const std::vector<std::vector<double>>&,
                                            const std::vector<std::vector<double>>&);

//Matrix specials
//Identity matrix
std::vector<std::vector<double>> mat_I(size_t);
//transpose matrix
std::vector<std::vector<double>> mat_T(const std::vector<std::vector<double>>&);
//companion matrix
std::vector<std::vector<double>> companion(const std::vector<double>&);

}
#endif // MATALG_H
