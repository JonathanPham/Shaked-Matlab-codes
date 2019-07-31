function [x,niter]=conj_multidiag_pre(A,x,b,nmax,tol)
%note that b and x must be column vectors
tic
lfil=70;
M=spalloc(nmax,nmax,lfil*nmax);
for j=1:nmax
     M(:,j)=r_sparse_inverse(A,j,tol,nmax);
     %M(:,j)=new_r_sparse_inverse(A,j,tol,nmax,lfil);
end
M=(M+M')/2;
toc;
tic;
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
toc;
end