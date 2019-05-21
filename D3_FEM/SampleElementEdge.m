function [N,h,je] = SampleElementEdge(s,r,e,k)
global IEN  Coord  nDoF  nNodesElement  nEdgesElement; 
%edgeIEN =; create some clever mapping
len=length(edgeIEN);
x = Coord(IEN(edgeIEN,e),1);
y = Coord(IEN(edgeIEN,e),2);
z = Coord(IEN(edgeIEN,e),3);
[Ne dNe]=buildnlin_bound(s,r);
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
% je= some function of dNe,x,y,z
end