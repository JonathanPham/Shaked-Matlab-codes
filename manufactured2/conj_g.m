function [x,niter]=conj_g(A,x,b,nmax,tol) %checked this function, it works!
%note that b and x must be column vectors
resn=b-A*x;
p=resn;
for niter=1:nmax
    res=resn;
    alpha=res'*res/(p'*A*p);
    x=x+alpha*p;
    resn=res-alpha*A*p;
    if norm(resn)<tol
        break;
    end    
    beta=(resn'*resn)/(res'*res);
    p=resn+beta*p;
end
end