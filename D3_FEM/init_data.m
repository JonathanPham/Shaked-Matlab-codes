function init_data(elementtype)
global nNodes nElements nNodesElement nDoF nEdgesElement ...
    Coord EBC NBC IEN f g h C Params grav faces;
%parameters
tol=1e-8;
L  = Params.L;
c  = Params.c;
t=   Params.t;
E = Params.E; 
v = Params.v;
grav=-1;
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
if (strcmpi(elementtype,'hex'))
    X=X(IEN(1:8,:));
    Y=Y(IEN(1:8,:));
    Z=Z(IEN(1:8,:));
    flag=4;
end
if (strcmpi(elementtype,'tet'))
    X=X(IEN(1:4,:));
    Y=Y(IEN(1:4,:));
    Z=Z(IEN(1:4,:));
    flag=3;
end
%find 6 faces of the body
H1 = abs(X)<tol;
H2 = abs(L-X) <tol;
H3 = abs(Y+c)<tol;
H4 = abs(Y-c)<tol;
H5 = abs(Z+t)<tol;
H6 = abs(Z-t)<tol;
% % Faces on NBC 
for e=1:nElements
    for i=1:nEdgesElement
        fvec=faces(1:flag,i);
        if (H5(fvec,e)==1)
            NBC(e,i)=5;
        elseif (H6(fvec,e)==1)
            NBC(e,i)=6;
        end
    end
end

C(:,1) = E;
C(:,2) = v;
end