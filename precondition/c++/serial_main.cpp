#include "matvecops.hpp"
#include "inv_sp.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

int main()
{
	std::vector<double> value, MM, m;
	std::vector<int> col, row, Mcol, Mrow;
	int k,j,count=0;
	const int nn=16; //size of Laplacian matrix
	const double tol=1.0e-8; // tolerance 
	// Define A in terms of CSC - this part has been verified
	value.push_back(2.0);
	col.push_back(0);
	row.push_back(0);
	for (j=0;j<nn-1;j++)
	{
		col.push_back(3*j+2);
		value.push_back(-1.0);
		row.push_back(j+1);
		value.push_back(-1.0);
		row.push_back(j);
		value.push_back(2.0);
		row.push_back(j+1);
	}
	col.push_back(3*nn-2);
	// solve m_j*a=e_j sparsely
	// MM.resize(n+ndiag*());
	Mcol.push_back(0);
	for (j=0;j<nn;j++)
	{
		m=inv_sp(value,col,row,j,nn,tol); //calculate approximate inverse's current row
		for (k=ndiag;k>=-ndiag;k--)
		{
			if (j-k>=0 && j-k<nn)
			{
				MM.push_back(m[j-k]); //save results of M
				Mrow.push_back(j-k); // save M's row indices
				count++;
			}
		} 
		Mcol.push_back(count);
		// MM[3*j]=m[j];
		// if (j>0)
		// 	MM[3*j-1]=m[j-1];
		// if (j<nn-1)
		// 	MM[3*j+1]=m[j+1];
	}
	// print sanity check
	for (j=0;j<Mcol.size();j++)
		std::cout<<Mcol[j]<<std::endl;
}