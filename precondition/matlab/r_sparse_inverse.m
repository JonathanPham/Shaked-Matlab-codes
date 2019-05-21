function [m]=r_sparse_inverse(A,j,tol,n)
    ej=spalloc(n,1,1);
    ej(j)=1;
    m=ej/A(j,j);
    r=ej-A*m;
    for k=1:20*n
        %d=sparsify_tol(r,10^-2,n,A);
        d=sparsify_shape(r,j,n,1);
        d(j)=r(j);
        %d=r;
        q=A*d;
        alpha=(r'*q)/(q'*q);
        m=m+alpha*d;
        r=r-alpha*q;
        if norm(r)<tol
            break;
        end
    end
end