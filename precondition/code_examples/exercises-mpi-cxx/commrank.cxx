/****************************************************************
 ****
 **** This program file is part of the book 
 **** `Parallel programming with MPI and OpenMP'
 **** by Victor Eijkhout, eijkhout@tacc.utexas.edu
 ****
 **** copyright Victor Eijkhout 2012-7
 ****
 **** MPI Exercise for the use of Comm_rank/size
 ****
 ****************************************************************/

#include <iostream>
#include <sstream>
using namespace std;
#include <mpi.h>

int main() {
  MPI_Comm comm = MPI_COMM_WORLD;
  int nprocs, procno, junk,len;
char name[50];
  
  MPI_Init(0,0);

  // Exercise part 1:
  // -- Use the routine
  //    MPI_Comm_size and MPI_Comm_rank
  // -- Let each processor print out a message like
  //    "Hello from processor 7 out of 12"
  //    reporting its number and the total number.
  stringstream proctext;
/**** your code here ****/
junk=MPI_Comm_size(comm,&nprocs);
junk=MPI_Comm_rank(comm,&procno);
junk=MPI_Get_processor_name(name,&len);
cout <<"Hello from processor "<< procno << " out of "<< nprocs << " running on "<< name<< endl;
  cout << proctext.str();

  // Exercise part 2:
  // -- let only processs zero print out
  //    "There are 16 processes"
  //    reporting only the total number
/**** your code here ****/
if(procno==0)
{
cout<< "There are " << nprocs <<" processes"<<endl;
}  
  MPI_Finalize();
  return 0;
}
