%% gmres
clear all;
clc;
n=2^20;
A=gallery('tridiag',(1:n-1).^2/2,1+(1:n).^2,(1:n-1).^2/2);
b=sparse(1:n,1,1);
tol=eps(1);
%tol=10^-8;
offdiag=2;
%[y,niter_y]=conj_g(A,b,b,n,tol);
%[x,niter]=conj_diag_pre(A,b,b,n,tol);
%tic;
[x,niter]=conj_multidiag_pre(A,b,b,n,tol,offdiag);
%toc;
%% minres
clear all;
n=10^4;
b=ones(n,1);
tol=eps(1);
A=gallery('tridiag',(1:n-1).^2/2,1+[1:n],(1:n-1).^2/2);
offdiag=0;
M=spalloc(n,n,(offdiag*2+1)*n);
for j=1:n
    M(:,j)=r_sparse_inverse(A,j,tol,n,offdiag);
end
M=(M+M')/2;
[ x, istop, itn, rnorm, Arnorm, Anorm, Acond, ynorm ] =  minres( A, b, M, 0, 0, 0, n, tol );
%% sparse inverse timing
clear all;
n=2^20;
tol=eps(1);
A=gallery('tridiag',(1:n-1).^2/2,1+[1:n],(1:n-1).^2/2);
offdiag=1;
M=spalloc(n,n,(offdiag*2+1)*n);
j=5;
M(:,j)=r_sparse_inverse(A,j,tol,n,offdiag);