#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

#include <vector>
/* dot product of two vectors */
double dotprod(std::vector<double> x,std::vector<double> y);

// vector addition/subtraction
std::vector<double> vecadd(std::vector<double> x,std::vector<double> y
    ,double alpha);

// scalar vector multiplication
std::vector<double> scalmult(double alpha,std::vector<double> x);


// Matrix (in CSC format) vector product
std::vector<double> matmult(std::vector<double> value,
    std::vector<int> row,std::vector<int> col, std::vector<double> x);

#endif /* MATVECOPS_HPP */
