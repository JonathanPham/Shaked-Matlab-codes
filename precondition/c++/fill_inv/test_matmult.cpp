#include "matvecops.hpp"
#include "inv_sp.hpp"
#include <vector>
#include <cmath>

int main()
{
	std::vector<double> value, M, m;
	std::vector<int> col, row;
	int j;
	const int nn=5; //size of Laplacian matrix
	// Define A in terms of CSC
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
	//
	m.push_back(-1);
	m.push_back(0);
	m.push_back(1);
	m.push_back(0.5);
	m.push_back(-0.5);
	m=matmult(value,col,row,m);
	for (j=0;j<nn;j++)
		std::cout<<m[j]<<std::endl;
}