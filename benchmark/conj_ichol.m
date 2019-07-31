function [x,niter]=conj_ichol(A,x,b,nmax,tol)
%note that b and x must be column vectors
tic
L=ichol(A);
nnz(L)/nmax
toc;
tic;
resn=b-A*x;
y=L\resn;
zn=L'\y;
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
    y=L\resn;
    zn=L'\y;
    beta=(zn'*resn)/(z'*res);
    p=zn+beta*p;
end
toc;
end