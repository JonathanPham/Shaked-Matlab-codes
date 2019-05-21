#ifndef _matvecops_HPP_
#define _matvecops_HPP_

#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>

//all functions here have been verified

/* dot product of two vectors */
double dotprod(std::vector<double> x,std::vector<double> y)
{
  double result=0;
  if (x.size()!=y.size())
  {
    std::cerr<<"ERROR: vectors must be of same size"<< std::endl;
  }
  for (unsigned i=0; i<x.size(); i++)
  {
    result+=x[i]*y[i];
  }
  return result;
}

// vector addition/subtraction (x+alpha*y) 
std::vector<double> vecadd(std::vector<double> x,std::vector<double> y
    ,double alpha)
{
  std::vector<double> result;
  if (x.size()!=y.size())
  {
    std::cerr<<"ERROR: vectors must be of same size"<< std::endl;
  }
  for (unsigned i=0; i<x.size(); i++)
  {
    result.push_back(x[i]+alpha*y[i]);
  }
  return result;
}
/* scalar vector multiplication */
std::vector<double> scalmult(double alpha,std::vector<double> x)
{
  std::vector<double> z;
  for (unsigned i=0; i<x.size(); i++)
  {
    z.push_back(x[i]*alpha);
  }
  return z;
}

/* Vector-Matrix (in CSC format) product */
std::vector<double> matmult(std::vector<double> value,
  std::vector<int> col,std::vector<int> row, std::vector<double> x)
{
  std::vector<double> result;
  for (unsigned i=0; i<col.size()-1;i++)
  {
    result.push_back(0);
    for (unsigned j=col[i];j<unsigned(col[i+1]);j++)
    {
      result[i]+=double(value[j]*x[int(row[j])]);
    }
    /* std::cout<<result[i]<<std::endl; */
  }
  return result;
}
#endif // _matvecops_HPP_