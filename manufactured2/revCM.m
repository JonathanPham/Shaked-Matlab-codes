function [R,R2]=revCM(A) %R2 is the inverse mapping to map back the solution
ind=A~=0;
A(ind)=1;
rsum=sum(A);
n=length(A);
%R2=zeros(1,n);
[~, ind_mind]=min(rsum);
R=zeros(1,n);
R2=zeros(1,n);
R(n)=ind_mind;
R2(ind_mind)=n;
i=1;
j=1;
while j<n
    tempR=find(A(R(n+1-i),:));
    adjR=setdiff(tempR,R);
    m=length(adjR);
    [~,sortR]=sort(rsum(adjR));
    R((n+1)-[j+1:j+m])=adjR(sortR);
    R2(adjR(sortR))=(n+1)-[j+1:j+m];
    i=i+1;
    j=j+m;
end
end