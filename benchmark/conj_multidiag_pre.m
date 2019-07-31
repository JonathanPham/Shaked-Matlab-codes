%% load matrix
clear all;
%load('bcsstk13.mat');
%load('crystm03.mat');
%load('aft01.mat');
load('finan512.mat');
%% problem setup
A=Problem.A;
tol=10^-8;
nmax=length(A);
b=sparse(A*ones(nmax,1));
x=spalloc(nmax,1,nmax);
%% sparse inverse
lfil=ceil(nnz(A)/nmax);
tic
M=entire_r_sparse_inverse(A,nmax,lfil);
M=(M+M')/2;
toc;
%% 
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