function [K,F]=diag_pre(A,x,b,nmax,tol)
vals=diag(K);
M=diag(1/sqrt(vals));
K=M*K*M;
F=M*K*M;
end