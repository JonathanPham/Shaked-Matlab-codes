function [d]=sparsify_tol(r,tol,n,A)
d=zeros(1,n);
for k=1:n
    rAk=A(k,:)*r';
    if((rAk^2)/(norm(A(k,:))^2)>tol)
        d(k)=r(k);
    end
end
end