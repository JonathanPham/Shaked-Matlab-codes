%% this file just has tests
L=2^7;
c=1;
t=1;
Nx=2;
Ny=1;
Nz=1;
x = linspace(0,L,Nx+1);
y = linspace(0,c,Ny+1);
z = linspace(0,t,Nz+1);
[X,Y,Z] = meshgrid(x,y,z);
X = permute(X,[2 1 3]);
Y = permute(Y,[2 1 3]);
Z=permute(Z,[2 1 3]);
XYZ=[X(:), Y(:), Z(:)];
TRIb = delaunayTriangulation(XYZ);
Coord=TRIb.Points;
IEN=(TRIb.ConnectivityList)';
%%
L=1;
c=1;
t=1;
Nx=2;
Ny=1;
Nz=1;
x = linspace(0,L,Nx+1);
y = linspace(0,c,Ny+1);
z = linspace(0,t,Nz+1);
[X,Y,Z] = meshgrid(x,y,z);
X = permute(X,[2 1 3]);
Y = permute(Y,[2 1 3]);
Z=permute(Z,[2 1 3]);
XYZ=[X(:), Y(:), Z(:)];
TRIb = delaunayTriangulation(XYZ);
Coord=TRIb.Points;
IEN=(TRIb.ConnectivityList)';
[Coord2, IEN2]=lin2quad(Coord, IEN);
%%
X=[1 0 0 0];
Y=[0 1 0 0];
Z=[0 0 1 0];
XYZ=[X(:), Y(:), Z(:)];
TRIb = delaunayTriangulation(XYZ);
Coord=TRIb.Points;
IEN=(TRIb.ConnectivityList)';
[Coord2, IEN2]=lin2quad(Coord, IEN);
%%
X=[1 0 0]';
Y=[0 1 0]';
TRIb = delaunayTriangulation(X,Y);
Coord=TRIb.Points;
IEN=(TRIb.ConnectivityList)';
[Coord2, IEN2]=lin2quad(Coord, IEN);
%% quad tet shape functions calculations
syms xi eta zeta alp Ne;
alp=1-xi-eta-zeta;
Ne=[xi*(2*xi-1) eta*(2*eta-1) zeta*(2*zeta-1) alp*(2*alp-1) 4*xi*eta 4*xi*zeta 4*xi*alp 4*eta*zeta 4*eta*alp 4*zeta*alp];
