function [x,niter]=conj_multidiag_pre(A,x,b,nmax,tol,offdiag)
%note that b and x must be column vectors
M=spalloc(nmax,nmax,(2*offdiag+1)*nmax);
for j=1:nmax
     M(j,:)=l_sparse_inverse(A,j,tol,nmax,offdiag);
end
M=(M+M')/2;
resn=b-A*x;
zn=M*resn;
p=zn;
for niter=1:nmax
    res=resn;
    z=zn;
    alpha=res'*z/(p'*A*p);
    x=x+alpha*p;
    resn=res-alpha*A*p;
    if norm(resn)<tol
        break;
    end
    zn=M*resn;
    beta=(zn'*resn)/(z'*res);
    p=zn+beta*p;
end
end