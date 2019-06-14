function [m]=l_sparse_inverse(A,j,tol,n,offdiag)
%ej=sparse(1,j,1,1,n,2*offdiag+1);
m=sparse(1,j,1/A(j,j),1,n,2*offdiag+1);
r=sparse(1,j,1,1,n,2*offdiag+1)-m*A;
for k=1:2*offdiag
    d=sparsify_shape(r,j,n,offdiag);
    q=d*A;
    if norm(q)<tol
        break;
    end
    alpha=(r*q')/(q*q');
    m=m+alpha*d;
    r=r-alpha*q;
    if norm(r)<tol
        break;
    end
end
end
