function init_data()
global nNodes nElements nNodesElement nDoF nEdgesElement ...
    Coord EBC NBC IEN f g h C Params grav;
%parameters
tol=1e-8;
L  = Params.L;
c  = Params.c;
t=   Params.t;
E = Params.E; 
v = Params.v;
grav=-9.8;
%allocation
C   = zeros(nElements,2); 
% f   = spalloc(nElements,nDoF,0);
f= [zeros(nElements,1) grav*ones(nElements,1) zeros(nElements,1)]/L/c/t/4;
g   = spalloc(nNodes,nDoF,round(nNodes/4)); 
EBC = spalloc(nNodes,nDoF,round(nNodes/4));
NBC = spalloc(nElements,nEdgesElement,round(nElements/4));  

% Step 1 set EBCs
X=Coord(:,1);
Y=Coord(:,2);
Z=Coord(:,3);

A=find(abs(X)<tol);
for i = 1:length(A)
    a = A(i);   
    EBC(a,1) = 1;
    EBC(a,2) = 1;
    EBC(a,3) =1; 
end

% Step 2 set NBCs - still needs work
%h = spalloc(nEdgesElement,nElements,round(nElements/4));
% X=X(IEN(1:8,:));
% %nodes on NBC
% H1 = abs(X)<tol;
% H2 = abs(L-X) <tol;
% % Faces on NBC - not sure how to define
% for e=1:nElements
%     for i=1:nEdgesElement
%         j=mod(i,nEdgesElement)+1; % finds the next face or goes to first if at last
%     end
% end

C(:,1) = E;
C(:,2) = v;
end