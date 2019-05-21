#include "matvecops.hpp"
#include "inv_sp.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <mpi.h>

int main()
{
	const int nn=2; //size of Laplacian matrix
	const double tol=1.0e-8; // tolerance 
	std::vector<double> value(3*nn-2,-1.0), m;
	std::vector<int> col(nn+1,0), row(3*nn-2,0), Mrow(nn+1,0);
	int k,count=0;
	int nprocs, procno;
	MPI_Comm comm = MPI_COMM_WORLD;
	
	
	// Define A in terms of CSC - 

	value[0]=2.0;

	for (k=0;k<nn-1;k++) //this loop might not need to be parallelized (relatively little time)
	{
		col[k+1]=3*k+2;
		//value[3*j+1]=-1.0;
		row[3*k+1]=k+1;
		//value[3*j+2]=-1.0;
		row[3*k+2]=k;
		value[3*k+3]=2.0;
		row[3*k+3]=k+1;
		if(k-ndiag>=0 && k+ndiag<nn) //middle of the matrix
			count+=2*ndiag+1;
		if(k-ndiag<0) //beginning of the matrix
			count+=1+ndiag+k;
		if(k+ndiag>=nn) // end of the matrix
			count+=ndiag+nn-k;
		Mrow[k+1]=count;
	}
	count+=1+ndiag;
	col[nn]=3*nn-2;
	Mrow[nn]=count;
	std::vector<int> Mcol(count,0);
	std::vector<double> MM(count,0);

	//solve m_j*a=e_j sparsely
	MPI_Init(NULL,NULL);
	MPI_Comm_size(comm,&nprocs);
	MPI_Comm_rank(comm,&procno);
	// nprocs=1;
	// procno=0;
	for (int j=procno;j<nn;j+=nprocs) //have to parallelize this
	{
		int i=0;
		m=inv_sp(value,col,row,j,nn,tol); //calculate approximate inverse's current row
		for (int p=ndiag;p>=-ndiag;p--)
		{
			if (j-p>=0 && j-p<nn)
			{
				//MM.push_back(m[j-k]); //save results of M
				MM[Mrow[j]+i]=m[j-p];
				//Mcol.push_back(j-k); // save M's col indices
				Mcol[Mrow[j]+i]=j-p;
				// std::cout<<"procno="<<procno<<std::endl;
				// std::cout<<"value="<<MM[Mrow[j]+i]<<std::endl;
				i++;
			}
		} 

	}
	// print sanity check
	// for (j=0;j<Mcol.size();j++)
	// 	std::cout<<Mcol[j]<<std::endl;
	if (procno==0)
		for (k=0;k<MM.size();k++)
		{
			std::cout<<MM[k]<<std::endl;
			std::cout<<Mcol[k]<<std::endl;
		}
	MPI_Finalize();
	return 0;
}