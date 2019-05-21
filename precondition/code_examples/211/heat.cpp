#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "heat.hpp"
#include "sparse.hpp"
#include "CGSolver.hpp"


int HeatEquation2D:: Setup(std::string inputfile)
{
  int n;
  double len, wid, h, tc, th;
  int i,j;//row, column index
  std::ifstream infile(inputfile.c_str());
  if (infile.is_open())
  {
    while (infile >> len >> wid >> h >> tc >> th);
  }
  else
  {
    std::cout << "Invalid file name"<<std::endl;
    std::terminate();
  }
  nx=(int) (len/h);//columns of solution
  ny=(int) (wid/h-1);//rows of solution
  for(i=0;i<nx;i++)
  {
    highb.push_back(th);
    /* lowb.push_back(tc*(2-exp(-10*pow(i*h-len/2,2)))); */
  }
  n=nx*ny;
  A.Resize(n,n);
  for(i=0;i<n;i++)
  {
    x.push_back(0.0);
    b.push_back(0.0);
    A.AddEntry(i,i,4);//diagonal points
  }
  for (i=0;i<ny;i++)//down points
  {
    A.AddEntry(i,ny*(nx-1)+i,-1);
    for(j=1;j<nx;j++)
    {
      A.AddEntry(j*ny+i,(j-1)*ny+i,-1);
    }
  }
  for(i=0;i<ny;i++)//right points
  {
    for (j=0;j<nx-1;j++)
    {
      A.AddEntry(j*ny+i,(j+1)*ny+i,-1);
    }
    A.AddEntry(ny*(nx-1)+i,i,-1);
  }
  for (j=0;j<nx;j++)//up points
  {
    b[j*ny]=th;
    for(i=1;i<ny;i++)
    {
      A.AddEntry(j*ny+i,j*ny+i-1,-1);
    }
  }
  //down points
  double cold;
  for(j=0;j<nx;j++)
  {
    for(i=0;i<ny-1;i++)
    {
      A.AddEntry(j*ny+i,j*ny+i+1,-1);
    }
    cold=tc*(2-exp(-10*pow(j*h-len/2,2)));
    lowb.push_back(cold);
    b[(j+1)*ny-1]=cold;
  }
  //Wrong code, ignore
  /* for(idxr=0;idxr<sizeA;idxr++) */
  /* { */
  /*   x.push_back(0); */
  /*   b.push_back(0); */
  /*   A.AddEntry(idxr,idxr,4); */
  /*   if(idxr%col==0) //we are on a left mesh point */
  /*   { */
  /*     A.AddEntry(idxr,idxr+col-1,-1); */
  /*   } */
  /*   else */
  /*   { */
  /*     A.AddEntry(idxr,idxr-1,-1); */
  /*   } */
  /*   if(idxr%col==col-1)// we are on a right mesh point */
  /*   { */
  /*     A.AddEntry(idxr,idxr-col+1,-1); */
  /*   } */
  /*   else */
  /*   { */
  /*     A.AddEntry(idxr,idxr+1,-1); */
  /*   } */
  /*   if(idxr<col) // we are in the first row bcs apply */
  /*   { */
  /*     b[idxr]+=th; */
  /*   } */
  /*   else */
  /*   { */
  /*     A.AddEntry(idxr,idxr-col,-1); */
  /*   } */
  /*   if(idxr+col>=sizeA)// we are in the last row bcs apply */
  /*   { */
  /*     b[idxr]+=tc*(2-exp(-10*pow((sizeA-idxr)*h-wid/2,2))); */
  /*   } */
  /*   else */
  /*   { */
  /*     A.AddEntry(idxr,idxr+col,-1); */
  /*   } */
  /* } */
  A.ConvertToCSR();
  return 0;
}

/* Method to solve system using CGsolver */
int HeatEquation2D::Solve(std::string soln_prefix)
{
  const double tol=0.00001;
  int niter;
  niter=CGSolver(A.get_a(),A.get_i_idx(),A.get_j_idx(),b,x,nx,ny,
      tol,soln_prefix,lowb,highb);
  /* std::cout<<"here"<<std::endl; */
  if (niter==-1) return niter;
  return 0;
}

/* TODO: Add any additional public methods you need */

