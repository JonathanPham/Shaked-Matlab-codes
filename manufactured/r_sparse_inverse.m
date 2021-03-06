function [m]=r_sparse_inverse(A,j,tol,n)
%ej=sparse(j,1,1,n,1,2*offdiag+1);
vecj=(abs(A(:,j))>tol);
sumj=sum(vecj);
m=sparse(j,1,1/A(j,j),n,1,sumj);
r=sparse(j,1,1,n,1)-A*m;
d=spalloc(n,1,sumj);
for k=1:sumj
    %d=r_sparsify_shape(r,j,n,1);
    d(vecj)=r(vecj);
    q=A*d;
    if norm(q)<tol
        break;
    end
    alpha=(r'*q)/(q'*q);
    m=m+alpha*d;
    r=r-alpha*q;
    if norm(r)<tol
        break;
    end
end
end