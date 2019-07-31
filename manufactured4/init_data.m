function init_data(elementtype)
global nNodes nElements nNodesElement nDoF nEdgesElement ...
    Coord EBC NBC IEN f g h C Params faces;
%parameters
tol=1e-8;
L  = Params.L;
c  = Params.c;
t=   Params.t;
E = Params.E; 
v = Params.v;
grav=Params.grav;
%allocation
C   = zeros(nElements,2); 
% f   = spalloc(nElements,nDoF,0);
%f= [zeros(nElements,1) grav*ones(nElements,1) zeros(nElements,1)]/L/c/t/4;
g   = spalloc(nNodes,nDoF,round(nNodes/4)); 
EBC = spalloc(nNodes,nDoF,round(nNodes/4));
NBC = spalloc(nElements,nEdgesElement,round(nElements/4));  

% Step 1 set EBCs
X=Coord(:,1);
Y=Coord(:,2);
Z=Coord(:,3);

A=find(abs(L-X)<tol);
A1=find(abs(c-Y)<tol);
A2=find(abs(t-Z)<tol);
A3=find(abs(L+X)<tol);
A4=find(abs(c+Y)<tol);
A5=find(abs(t+Z)<tol);
EBC(A,:)=1;
EBC(A1,:)=1;
EBC(A2,:)=1;
EBC(A3,:)=1;
EBC(A4,:)=1;
EBC(A5,:)=1;

C(:,1) = E;
C(:,2) = v;
end