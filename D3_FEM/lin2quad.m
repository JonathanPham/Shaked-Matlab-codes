function [Coord, IENnew]=lin2quad(Coord,IEN)
nNodes=size(Coord,1);
nElements=size(IEN,2);
newn=6;
ndim=3;
nNodesElement=size(IEN,1);
Coord = [ Coord; zeros(newn*nElements,ndim) ];
IENnew   = [ IEN;   zeros(newn,nElements)   ];
mID   = spalloc(nNodes,nNodes,newn*nElements);
% Step 2: Add midside nodes to all edges
for e = 1:nElements
    AB=nchoosek(IEN(:,e),2);
    for i = 1:newn
        
        % Step 2a: Get global number for nodes defining the edge
        A = AB(i,1);
        B = AB(i,2);
        
        % Step 2b: If edge doesn't have a node, add it
        if mID(A,B) == 0
            
            % Assign global node number
            nNodes     = nNodes + 1;
            mID(A,B)   = nNodes;
            mID(B,A)   = nNodes;
            IENnew(nNodesElement+i,e) = nNodes;
            
            xm = 1/2*(Coord(A,:) + Coord(B,:));
            Coord(nNodes,:) = xm; 
            % Step 2c: If node already on edge, update IEN
        else
            IENnew(i+nNodesElement,e) = mID(A,B);
        end
        
    end
end

% Step 3: Remove extra entries in Coord
Coord = Coord(1:nNodes,:);
end