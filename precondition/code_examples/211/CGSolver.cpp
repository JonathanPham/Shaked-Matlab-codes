#ifndef CGSOLVER_HPP
#define CGSOLVER_HPP
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "matvecops.hpp"
/* Function that implements the CG algorithm for a linear system
 *  *
 *   * Ax = b
 *    *
 *     * where A is in CSR format.  The starting guess for the solution
 *      * is provided in x, and the solver runs a maximum number of iterations
 *       * equal to the size of the linear system.  Function returns the
 *        * number of iterations to converge the solution to the specified
 *         * tolerance, or -1 if the solver did not converge.
 *          */

void printx(std::vector<double> &x,std::vector<double> &lowb,
    std::vector<double> &highb, int niter,int nx, int ny,
    std::string soln_prefix)
{
  std::string outstr=soln_prefix;
  int i,j;
  if (niter<10) outstr+="0";
  if (niter<100) outstr+="0";
  outstr+=std::to_string(niter);
  outstr+=".txt";
  std::ofstream outfile(outstr.c_str());
  if (outfile.is_open())
  {
    std::setprecision(4);
    for (i=0 ;i<nx;i++)//top row
    {
      outfile<<highb[i]<<"\n"<<std::scientific;
    }
    for(i=0;i<ny;i++)//x solution
    {
      for(j=0;j<nx;j++)
      {
        outfile<<x[j*ny+i]<<"\n"<<std::scientific;
      }
    }
    for (i=0;i<nx;i++)// bottom row
    {
      outfile<<lowb[i]<<"\n"<<std::scientific;
    }
    outfile.close();
  }
  else
  {
    std::cerr<<"ERROR: failed to open file"<< std::endl;
  }
}
int CGSolver(std::vector<double> &val,
    std::vector<int>    &row_ptr,
    std::vector<int>    &col_idx,
    std::vector<double> &b,
    std::vector<double> &x,int nx, int ny,
    double              tol,
    std::string soln_prefix,
    std::vector<double> &lowb,
    std::vector<double> &highb)
{
  std::vector<double> res;
  std::vector<double> p;
  double alpha,bet,normc,norml,norm0;/*I will be working with squared norms*/
  int niter=0,nitermax;
  nitermax=int(row_ptr.size());
  res=vecadd(b,matmult(val,row_ptr,col_idx,x),-1);
  norm0=dotprod(res,res);/*this is actually norm squared*/
  p=res;/*this is intentionally copying the residual to p*/
  while(niter<nitermax)
  {
    alpha=dotprod(res,res)/dotprod(p,matmult(val,row_ptr,col_idx,p));
    x=vecadd(x,scalmult(alpha,p),1);
    niter++;
    norml=dotprod(res,res);/*I only care about the norm of the residual*/
    res=vecadd(res,scalmult(alpha,matmult(val,row_ptr,col_idx,p)),-1);
    normc=dotprod(res,res);
    if (normc/norm0<tol*tol)
    {
      printx(x,lowb,highb,niter-1,nx,ny,soln_prefix);
      std::cout<<"SUCCESS: CG solver converged in "<<niter << " iterations."
        <<std::endl;
      return niter;
    }
    if (niter%10==1)
    {
      printx(x,lowb,highb,niter-1,nx,ny,soln_prefix);
    }
    bet=normc/norml;
    p=vecadd(res,scalmult(bet,p),1);
  }
  return -1;
}

#endif /* CGSOLVER_HPP */
