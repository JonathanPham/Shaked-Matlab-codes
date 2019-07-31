#ifndef _INV_SP_HPP_
#define _INV_SP_HPP_

#define ndiag 1 //number of off diagonals per side

#include "matvecops.hpp"
#include <vector>
#include <cmath>

//correct!

std::vector<double> inv_sp(std::vector<double> &val,
    std::vector<int>    &col_ptr,
    std::vector<int>    &row_idx,int j, int n, double tol)
{
	int i,k;
	double alpha;
	std::vector<double> ej(n,0),m(n,0),r,q,d(n,0);
	// ej.resize(n,0);
	// m.resize(n,0);
	m[j]=1.0/val[3*j];
	// std::cout<<m[0]<< " " << m[1]<<std::endl;
	// d.resize(n,0);
	ej[j]=1.0;
	r=vecadd(ej,matmult(val,col_ptr,row_idx,m),-1.0);
	// std::cout<<r[0]<< " " << r[1]<<std::endl;
	for (k=0;k<=ndiag*2;k++) //gives us as many DOF as non-zero entries
	{
		//This part of the code sparsifies d (effectively assumes banded structure)
		d[j]=r[j];
		for (i=1;i<=ndiag;i++)
		{
			if (j>=i)
				d[j-i]=r[j-i];
			if (j<n-i)
				d[j+i]=r[j+i];
		}
		//
		// std::cout<<d[0]<< " " << d[1]<<std::endl;
		q=matmult(val,col_ptr,row_idx,d);
		// std::cout<<q[0]<< " " << q[1]<<std::endl;
		alpha=dotprod(r,q)/dotprod(q,q);
		// std::cout<<alpha<<std::endl;
		//m+=scalmult(alpha,d);
		m=vecadd(m,d,alpha);
		//r-=scalmult(alpha,q);
		r=vecadd(r,q,-alpha);
		if (dotprod(r,r)<tol)
			break;
	}
	return m;
}

#endif // _INV_SP_HPP_
