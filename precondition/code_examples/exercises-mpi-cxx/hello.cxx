/****************************************************************
 ****
 **** This program file is part of the book 
 **** `Parallel programming with MPI and OpenMP'
 **** by Victor Eijkhout, eijkhout@tacc.utexas.edu
 ****
 **** copyright Victor Eijkhout 2012-6
 ****
 **** MPI Exercise
 ****
 ****************************************************************/

#include <iostream>
#include <mpi.h>
#include "mpi.h"
using namespace std;

int main(int argc,char **argv) {
int ierr;

cout << "hello world!" << endl;
ierr=MPI_Init(&argc,&argv);

cout << "hello world!" << endl;
ierr=MPI_Finalize();
 
cout << "hello world!" << endl;
  return 0;
}
