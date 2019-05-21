% test for sparse inverse
n=2^15;
%n=8;
A = spalloc(n,n,3*n);
M = spalloc(n,n,11*n); %the inverse will be stored here
for i=1:n %this is just a test for a very simple A
    A(i,i)=2;
    if i>1
        A(i,i-1)=-1;
    end
    if i<n
        A(i,i+1)=-1;
    end
end
for j=1:n
    M(j,:)=l_sparse_inverse(A,j,10^-8,n);
end
% test this vs diagonal preconditioner
my_err=norm(speye(n)-M*A,'fro');
%%
W=spalloc(n,n,n);
for i=1:n
    W(i,i)=1/sqrt(A(i,i));
end
diag_err=norm(eye(n,n)-W*A*W,'fro');
%%
n=32768;
value=zeros(3*n-2,1);
col=zeros(n+1,1);
row=zeros(3*n-2,1);
value(1)=2;
row(1)=1;
col(1)=1;
for j=0:n-2
    col(j+2)=3*j+3;
    row(3*j+2)=j+2;
    value(3*j+2)=-1;
    row(3*j+3)=j+1;
    value(3*j+3)=-1;
    row(3*j+4)=j+2;
    value(3*j+4)=2;
end
col(n+1)=3*n-1;