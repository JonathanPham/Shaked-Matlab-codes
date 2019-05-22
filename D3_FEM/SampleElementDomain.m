function [N,B,je] = SampleElementDomain(e1,e2,e3,e,elementtype) % this works
global IEN  Coord  nNodesElement  nDoF;
X = Coord(IEN(:,e),:);
[Ne,dNe]=buildnlin_der(e1,e2,e3);
% if (e==1)
%     X
%     dNe
% end
J=dNe*X;
dNdx=J\dNe;
N=zeros(3,nDoF*nNodesElement);
B=zeros(6,nDoF*nNodesElement);
for j=1:nNodesElement
    p=(j-1)*nDoF+1;
    q=p+1;
    r=q+1;
    N(1,p)=Ne(j);
    N(2,q)=Ne(j);
    N(3,q)=Ne(j);
    B(1,p)=dNdx(1,j);
    B(2,q)=dNdx(2,j);
    B(3,r)=dNdx(3,j);
    B(4,q)=dNdx(3,j);
    B(4,r)=dNdx(2,j);
    B(5,p)=dNdx(3,j);
    B(5,r)=dNdx(1,j);
    B(6,p)=dNdx(2,j);
    B(6,q)=dNdx(1,j);
end
je=det(J);
end