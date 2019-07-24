load('bcsstk32.mat'); %fails to converge
%load('qa8fk.mat');
%load('onetone2.mat');
%load('epb3.mat');
%load('Wordnet3.mat');
%load('ex5.mat');
%load('bcsstm13.mat');
%load('olm500.mat');
A=Problem.A;
n=length(A);
b=sparse(1:n,1,1);
tol=10^-8;
%x=A\b;
%[x,niter]=conj_multidiag_pre(A,b,b,n,tol);
%[x,niter]=conj_ichol(A,b,b,n,tol);
M=spalloc(n,n,5*n);
for j=1:n
    M(j,:)=r_sparse_inverse(A,j,tol,n);
end
M=(M+M')/2;
[ d1, istop, itn, rnorm, Arnorm, Anorm, Acond, ynorm ] =  minres( A, b, M, 0, 0, 0, n, tol );
norm(A*d1-b);