function [N,h,je] = SampleElementEdge(s,r,e,k,elementtype,order)
global IEN  Coord  nDoF  nNodesElement  nEdgesElement faces; 
edgeIEN =faces(:,k);
len=length(edgeIEN);
x = Coord(IEN(edgeIEN,e),1);
y = Coord(IEN(edgeIEN,e),2);
z = Coord(IEN(edgeIEN,e),3);
[Ne, dNe]=buildnlin_bound(s,r,elementtype,order);
h=SampleEdgeTraction(Ne*x,Ne*y,Ne*z,e,k);
N=zeros(nDoF,nDoF*nNodesElement); %something here is funky about N
for j=1:len %bounds for j not clear with definitions of N, edgeIEN
    p=(j-1)*nDoF+1;
    q=p+1;
    r=q+1;
    N(1,p)=Ne(j);
    N(2,q)=Ne(j);
    N(3,r)=Ne(j);
end
xs=dNe(1,:)*x;
xr=dNe(2,:)*x;
ys=dNe(1,:)*y;
yr=dNe(2,:)*y;
zs=dNe(1,:)*z;
zr=dNe(2,:)*z;
je= norm(det([ys zs; yr zr])-det([xs zs; xr zr])+det([xs ys; xr yr]));
end