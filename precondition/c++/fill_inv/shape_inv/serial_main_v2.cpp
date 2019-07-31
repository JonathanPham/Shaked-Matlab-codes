#include "matvecops.hpp"
#include "inv_sp.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
// #include <mpi.h>

int main()
{
	const int nn=8; //size of Laplacian matrix
	const double tol=1.0e-8; // tolerance 
	std::vector<double> value(3*nn-2,-1.0), m;
	std::vector<int> col(nn+1,0), row(3*nn-2,0), Mrow(nn+1,0);
	int k,j,count=0,i;
	// Define A in terms of CSC - 

	value[0]=2.0;

	for (j=0;j<nn-1;j++) //this loop might not need to be parallelized (relatively little time)
	{
		col[j+1]=3*j+2;
		//value[3*j+1]=-1.0;
		row[3*j+1]=j+1;
		//value[3*j+2]=-1.0;
		row[3*j+2]=j;
		value[3*j+3]=2.0;
		row[3*j+3]=j+1;
		if(j-ndiag>=0 && j+ndiag<nn) //middle of the matrix
			count+=2*ndiag+1;
		if(j-ndiag<0) //beginning of the matrix
			count+=1+ndiag+j;
		if(j+ndiag>=nn) // end of the matrix
			count+=ndiag+nn-j;
		Mrow[j+1]=count;
	}
	count+=1+ndiag;
	std::cout<<count<<std::endl;
	col[nn]=3*nn-2;
	Mrow[nn]=count;
	std::vector<int> Mcol(count,0);
	std::vector<double> MM(count,0);

	// solve m_j*a=e_j sparsely

	for (j=0;j<nn;j++) //have to parallelize this
	{
		i=0;
		m=inv_sp(value,col,row,j,nn,tol); //calculate approximate inverse's current row
		for (k=ndiag;k>=-ndiag;k--)
		{
			if (j-k>=0 && j-k<nn)
			{
				//MM.push_back(m[j-k]); //save results of M
				MM[Mrow[j]+i]=m[j-k];
				//Mcol.push_back(j-k); // save M's col indices
				Mcol[Mrow[j]+i]=j-k;
				i++;
			}
		} 

		// MM[3*j]=m[j];
		// if (j>0)
		// 	MM[3*j-1]=m[j-1];
		// if (j<nn-1)
		// 	MM[3*j+1]=m[j+1];
	}

	// print sanity check
	// for (j=0;j<Mcol.size();j++)
	// 	std::cout<<Mcol[j]<<std::endl;
	for (j=0;j<MM.size();j++)
		std::cout<<MM[j]<<std::endl;
}