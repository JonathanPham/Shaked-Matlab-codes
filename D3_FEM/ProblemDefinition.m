function ProblemDefinition(elementtype,order) 
global nDim nNodes nElements nNodesElement nDoF nEquations ...
    nEdgesElement Coord ID IEN LM EBC Params f;
%Solution paramaters
nDoF=3;
nDim=3;
nEdgesElement=6; %this is actually the number of faces in 3D

%Mesh parameters
Params.Nx = 2^7;  % Number of elements along x-axis
Params.Ny = 2^1;  % Number of elements along y-axis
Params.Nz = 2^1;  % Number of elements along z-axis

% Parameters for Plate Geometry
Params.L = 2^7;      % Length of plate, i.e. (0 < x < L)
Params.c =  1;      % Half-height of plate, i.e. (-c < y < c)
Params.t =   1;      % Thickness of plate (0<z<t)

% Parameters for Material Properties
Params.E = 10^7;     % Youngs Modulus
Params.v = 0.25;     % Poisson Ratio

% Parameters for Loading Conditions
Params.P = 80;       % Shear load at beam end

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
if (strcmpi(elementtype,'hex'))
    Coord = [X(:), Y(:), Z(:)];
    % Build IEN Array
    IEN=quadIEN(Nx,Ny,Nz,order);
elseif (strcmpi(elementtype,'tet'))
    XYZ=[X(:), Y(:), Z(:)];
    TRI = delaunayTriangulation(XYZ);
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

init_data();

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