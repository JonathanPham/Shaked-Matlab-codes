function [d]=r_sparsify_shape(r,j,n,k)
d=spalloc(n,1,2*k+1);
d(j)=r(j);
for i=1:k
if j>i
    d(j-i)=r(j-i);
end
if j<=n-i
    d(j+i)=r(j+i);
end
end