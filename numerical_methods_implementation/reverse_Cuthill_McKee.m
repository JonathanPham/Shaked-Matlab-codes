%% first create a sparse symmetric matrix with random connectivity
n=2^7;
A = randi([0 1], n,n);
AA=A'*A;
AA=(~mod(AA,n/8));
AA=AA-diag(diag(AA));
% r=sum(AA);
% threshmin=0.65*min(diag(AA)); %chose a number that would maintain connectivity
% ind=(AA<threshmin);
% threshmax=0.7*min(diag(AA));
% ind2=(AA>threshmax);
% AA(ind)=0;
% AA(~ind)=1;
% AA(ind2)=0;
%sum(sum(AA))
%%
%AA=AA-diag(diag(AA));
%spAA=sparse(AA);
rsum=sum(AA);
[mind, ind_mind]=min(rsum);
R=ind_mind;
i=1;
while length(R)<n
    tempR=find(AA(R(i),:));
    adjR=setdiff(tempR,R);
%     if isempty(adjR)
%         break;
%     end
    [junk,sortR]=sort(rsum(adjR));
    R=[R adjR(sortR)];
    i=i+1;
end
R=flip(R);