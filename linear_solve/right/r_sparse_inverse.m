function [m]=r_sparse_inverse(A,j,tol,n,offdiag)
%ej=sparse(j,1,1,n,1,2*offdiag+1);
m=sparse(j,1,1/A(j,j),n,1,2*offdiag+1);
r=sparse(j,1,1,n,1,2*offdiag+1)-A*m;
for k=1:2*offdiag
    d=r_sparsify_shape(r,j,n,1);
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