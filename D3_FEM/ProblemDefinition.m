function ProblemDefinition(elementtype,order) 
global nDim nNodes nElements nNodesElement nDoF nEquations ...
    nEdgesElement Coord ID IEN LM EBC Params f faces;
%Solution paramaters
nDoF=3;
nDim=3;
if (strcmpi(elementtype,'hex'))
    nEdgesElement=6; %this is actually the number of faces in 3D
    switch order
        case 1
            faces=[1,2,3,4;1,5,6,2;1,4,8,5;2,3,7,6;3,4,8,7;5,6,7,8]';
        case 2
            faces=[1,2,3,4,9,10,11,12,21;1,5,6,2,17,13,18,9,23;1,4,8,5,12,20,16,17,25;...
                2,3,7,6,10,19,14,18,26;3,4,8,7,11,20,15,19,24;5,6,7,8,13,14,15,16,22]';
    end
end
if (strcmpi(elementtype,'tet'))
    nEdgesElement=4; %this is actually the number of faces in 3D
    switch order
        case 1
            faces=[1,2,3;1,2,4;1,4,3;2,3,4,]';
        case 2
            faces=[1,2,3,5,8,6;1,2,4,5,9,7;1,4,3,7,10,6;2,3,4,8,10,9]';
    end
end

%Mesh parameters
Params.Nx = 2^10;  % Number of elements along x-axis
Params.Ny = 2^1;  % Number of elements along y-axis
Params.Nz = 2^1;  % Number of elements along z-axis

% Parameters for Plate Geometry
Params.L = 2^7;      % Length of plate, i.e. (0 < x < L)
Params.c =  1;      % Half-height of plate, i.e. (-c < y < c)
Params.t =   1;      % Thickness of plate (0<z<t)

% Parameters for Material Properties
Params.E = 10^7;     % Youngs Modulus
Params.v = 0;     % Poisson Ratio

% Parameters for forces
Params.P = 0;       % traction on z boundaries
Params.grav=-1;

%Mesh paramaters
Nx = Params.Nx;
Ny = Params.Ny;
Nz = Params.Nz;
L  = Params.L;
c  = Params.c;
t=   Params.t;

% Build Mesh

if (strcmpi(elementtype,'hex'))
    x = linspace(0,L,order*Nx+1);
    y = linspace(-c,c,order*Ny+1);
    z = linspace(-t,t,order*Nz+1);
end
if (strcmpi(elementtype,'tet'))
    x = linspace(0,L,Nx+1);
    y = linspace(-c,c,Ny+1);
    z = linspace(-t,t,Nz+1);
end
[X,Y,Z] = meshgrid(x,y,z);
X = permute(X,[2 1 3]);
Y = permute(Y,[2 1 3]);
Z=permute(Z,[2 1 3]);
XYZ=[X(:), Y(:), Z(:)];
if (strcmpi(elementtype,'hex'))
    Coord = [X(:), Y(:), Z(:)];
    % Build IEN Array
    IEN=quadIEN(Nx,Ny,Nz,order);
    %patch('Vertices',Coord(),'Faces',IEN,'FaceVertexCData',hsv(1),'FaceColor','none');
elseif (strcmpi(elementtype,'tet'))
    TRI = delaunayTriangulation(XYZ);
    figure
    tetramesh(TRI);
    Coord=TRI.Points;
    IEN=(TRI.ConnectivityList)';
    if order==2
        [Coord, IEN]=lin2quad(Coord, IEN,elementtype);
    end
else
    error('unsupported element type');
end
nNodes=size(Coord,1);
nElements=size(IEN,2);
nNodesElement=size(IEN,1);

init_data(elementtype);

%Allocate arrays
ID=zeros(nDoF,nNodes);
LM=zeros(nNodesElement*nDoF,nElements);

%Create ID
I=full(EBC'==0);
nEquations=sum(sum(I));
ID(I)=1:nEquations;
ID=ID';

%Create LM
P=ID(IEN,:)';
LM(:)=P(:);
end