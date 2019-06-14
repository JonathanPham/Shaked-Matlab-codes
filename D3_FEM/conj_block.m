function [x,niter]=conj_block(A,x,b,nmax,tol,blk)
%note that b and x must be column vectors
if mod(nmax,blk)~=0
    error('nmax%blk must be equal to 0');
end
M=spalloc(nmax,nmax,blk*nmax);
for j=1:nmax/blk
    vec=(j-1)*blk+1:blk*j;
     M(vec,vec)=inv(A(vec,vec));
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