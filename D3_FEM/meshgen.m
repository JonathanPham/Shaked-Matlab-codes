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
%%
L=1;
c=1;
t=1;
Nx=1;
Ny=1;
Nz=1;
order=2;
x = linspace(-L,L,order*Nx+1);
y = linspace(-c,c,order*Ny+1);
z = linspace(-t,t,order*Nz+1);
[X,Y,Z] = meshgrid(x,y,z);
X = permute(X,[2 1 3]);
Y = permute(Y,[2 1 3]);
Z=permute(Z,[2 1 3]);
XYZ=[X(:), Y(:), Z(:)];
%% quad hex shape functions
syms xi eta zeta L1 L2 L3 Ne;
L1=[1/2*xi*(xi-1) 1/2*xi*(xi+1) 1-xi^2];
L2=subs(L1,xi,eta);
L3=subs(L1,xi,zeta);
Ne=[L1(1)*L2(1)*L3(1) L1(2)*L2(1)*L3(1) L1(2)*L2(2)*L3(1) L1(1)*L2(2)*L3(1)...
    L1(1)*L2(1)*L3(2) L1(2)*L2(1)*L3(2) L1(2)*L2(2)*L3(2) L1(1)*L2(2)*L3(2)...
    L1(3)*L2(1)*L3(1) L1(2)*L2(3)*L3(1) L1(3)*L2(2)*L3(1) L1(1)*L2(3)*L3(1)...
    L1(3)*L2(1)*L3(2) L1(2)*L2(3)*L3(2) L1(3)*L2(2)*L3(2) L1(1)*L2(3)*L3(2)...
    L1(1)*L2(1)*L3(3) L1(2)*L2(1)*L3(3) L1(2)*L2(2)*L3(3) L1(1)*L2(2)*L3(3)...
    L1(3)*L2(3)*L3(1) L1(3)*L2(3)*L3(2) L1(3)*L2(1)*L3(3) L1(3)*L2(2)*L3(3)...
    L1(1)*L2(3)*L3(3) L1(2)*L2(3)*L3(3) L1(3)*L2(3)*L3(3)];
dNe1=[ (eta*xi*zeta*(eta - 1)*(zeta - 1))/8 + (eta*zeta*(eta - 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(eta - 1)*(zeta - 1))/8 + (eta*zeta*(eta - 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(eta + 1)*(zeta - 1))/8 + (eta*zeta*(eta + 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(eta + 1)*(zeta - 1))/8 + (eta*zeta*(eta + 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(eta - 1)*(zeta + 1))/8 + (eta*zeta*(eta - 1)*(xi - 1)*(zeta + 1))/8, (eta*xi*zeta*(eta - 1)*(zeta + 1))/8 + (eta*zeta*(eta - 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(eta + 1)*(zeta + 1))/8 + (eta*zeta*(eta + 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(eta + 1)*(zeta + 1))/8 + (eta*zeta*(eta + 1)*(xi - 1)*(zeta + 1))/8, -(eta*xi*zeta*(eta - 1)*(zeta - 1))/2, - (xi*zeta*(eta^2 - 1)*(zeta - 1))/4 - (zeta*(eta^2 - 1)*(xi + 1)*(zeta - 1))/4, -(eta*xi*zeta*(eta + 1)*(zeta - 1))/2, - (xi*zeta*(eta^2 - 1)*(zeta - 1))/4 - (zeta*(eta^2 - 1)*(xi - 1)*(zeta - 1))/4, -(eta*xi*zeta*(eta - 1)*(zeta + 1))/2, - (xi*zeta*(eta^2 - 1)*(zeta + 1))/4 - (zeta*(eta^2 - 1)*(xi + 1)*(zeta + 1))/4, -(eta*xi*zeta*(eta + 1)*(zeta + 1))/2, - (xi*zeta*(eta^2 - 1)*(zeta + 1))/4 - (zeta*(eta^2 - 1)*(xi - 1)*(zeta + 1))/4, - (eta*xi*(zeta^2 - 1)*(eta - 1))/4 - (eta*(zeta^2 - 1)*(eta - 1)*(xi - 1))/4, - (eta*xi*(zeta^2 - 1)*(eta - 1))/4 - (eta*(zeta^2 - 1)*(eta - 1)*(xi + 1))/4, - (eta*xi*(zeta^2 - 1)*(eta + 1))/4 - (eta*(zeta^2 - 1)*(eta + 1)*(xi + 1))/4, - (eta*xi*(zeta^2 - 1)*(eta + 1))/4 - (eta*(zeta^2 - 1)*(eta + 1)*(xi - 1))/4, xi*zeta*(eta^2 - 1)*(zeta - 1), xi*zeta*(eta^2 - 1)*(zeta + 1), eta*xi*(zeta^2 - 1)*(eta - 1), eta*xi*(zeta^2 - 1)*(eta + 1), (xi*(eta^2 - 1)*(zeta^2 - 1))/2 + ((eta^2 - 1)*(zeta^2 - 1)*(xi - 1))/2, (xi*(eta^2 - 1)*(zeta^2 - 1))/2 + ((eta^2 - 1)*(zeta^2 - 1)*(xi + 1))/2, -2*xi*(eta^2 - 1)*(zeta^2 - 1)];
            dNe2=[ (eta*xi*zeta*(xi - 1)*(zeta - 1))/8 + (xi*zeta*(eta - 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(xi + 1)*(zeta - 1))/8 + (xi*zeta*(eta - 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(xi + 1)*(zeta - 1))/8 + (xi*zeta*(eta + 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(xi - 1)*(zeta - 1))/8 + (xi*zeta*(eta + 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(xi - 1)*(zeta + 1))/8 + (xi*zeta*(eta - 1)*(xi - 1)*(zeta + 1))/8, (eta*xi*zeta*(xi + 1)*(zeta + 1))/8 + (xi*zeta*(eta - 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(xi + 1)*(zeta + 1))/8 + (xi*zeta*(eta + 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(xi - 1)*(zeta + 1))/8 + (xi*zeta*(eta + 1)*(xi - 1)*(zeta + 1))/8, - (eta*zeta*(xi^2 - 1)*(zeta - 1))/4 - (zeta*(xi^2 - 1)*(eta - 1)*(zeta - 1))/4, -(eta*xi*zeta*(xi + 1)*(zeta - 1))/2, - (eta*zeta*(xi^2 - 1)*(zeta - 1))/4 - (zeta*(xi^2 - 1)*(eta + 1)*(zeta - 1))/4, -(eta*xi*zeta*(xi - 1)*(zeta - 1))/2, - (eta*zeta*(xi^2 - 1)*(zeta + 1))/4 - (zeta*(xi^2 - 1)*(eta - 1)*(zeta + 1))/4, -(eta*xi*zeta*(xi + 1)*(zeta + 1))/2, - (eta*zeta*(xi^2 - 1)*(zeta + 1))/4 - (zeta*(xi^2 - 1)*(eta + 1)*(zeta + 1))/4, -(eta*xi*zeta*(xi - 1)*(zeta + 1))/2, - (eta*xi*(zeta^2 - 1)*(xi - 1))/4 - (xi*(zeta^2 - 1)*(eta - 1)*(xi - 1))/4, - (eta*xi*(zeta^2 - 1)*(xi + 1))/4 - (xi*(zeta^2 - 1)*(eta - 1)*(xi + 1))/4, - (eta*xi*(zeta^2 - 1)*(xi + 1))/4 - (xi*(zeta^2 - 1)*(eta + 1)*(xi + 1))/4, - (eta*xi*(zeta^2 - 1)*(xi - 1))/4 - (xi*(zeta^2 - 1)*(eta + 1)*(xi - 1))/4, eta*zeta*(xi^2 - 1)*(zeta - 1), eta*zeta*(xi^2 - 1)*(zeta + 1), (eta*(xi^2 - 1)*(zeta^2 - 1))/2 + ((xi^2 - 1)*(zeta^2 - 1)*(eta - 1))/2, (eta*(xi^2 - 1)*(zeta^2 - 1))/2 + ((xi^2 - 1)*(zeta^2 - 1)*(eta + 1))/2, eta*xi*(zeta^2 - 1)*(xi - 1), eta*xi*(zeta^2 - 1)*(xi + 1), -2*eta*(xi^2 - 1)*(zeta^2 - 1)];
            dNe3=[ (eta*xi*zeta*(eta - 1)*(xi - 1))/8 + (eta*xi*(eta - 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(eta - 1)*(xi + 1))/8 + (eta*xi*(eta - 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(eta + 1)*(xi + 1))/8 + (eta*xi*(eta + 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(eta + 1)*(xi - 1))/8 + (eta*xi*(eta + 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(eta - 1)*(xi - 1))/8 + (eta*xi*(eta - 1)*(xi - 1)*(zeta + 1))/8, (eta*xi*zeta*(eta - 1)*(xi + 1))/8 + (eta*xi*(eta - 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(eta + 1)*(xi + 1))/8 + (eta*xi*(eta + 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(eta + 1)*(xi - 1))/8 + (eta*xi*(eta + 1)*(xi - 1)*(zeta + 1))/8, - (eta*zeta*(xi^2 - 1)*(eta - 1))/4 - (eta*(xi^2 - 1)*(eta - 1)*(zeta - 1))/4, - (xi*zeta*(eta^2 - 1)*(xi + 1))/4 - (xi*(eta^2 - 1)*(xi + 1)*(zeta - 1))/4, - (eta*zeta*(xi^2 - 1)*(eta + 1))/4 - (eta*(xi^2 - 1)*(eta + 1)*(zeta - 1))/4, - (xi*zeta*(eta^2 - 1)*(xi - 1))/4 - (xi*(eta^2 - 1)*(xi - 1)*(zeta - 1))/4, - (eta*zeta*(xi^2 - 1)*(eta - 1))/4 - (eta*(xi^2 - 1)*(eta - 1)*(zeta + 1))/4, - (xi*zeta*(eta^2 - 1)*(xi + 1))/4 - (xi*(eta^2 - 1)*(xi + 1)*(zeta + 1))/4, - (eta*zeta*(xi^2 - 1)*(eta + 1))/4 - (eta*(xi^2 - 1)*(eta + 1)*(zeta + 1))/4, - (xi*zeta*(eta^2 - 1)*(xi - 1))/4 - (xi*(eta^2 - 1)*(xi - 1)*(zeta + 1))/4, -(eta*xi*zeta*(eta - 1)*(xi - 1))/2, -(eta*xi*zeta*(eta - 1)*(xi + 1))/2, -(eta*xi*zeta*(eta + 1)*(xi + 1))/2, -(eta*xi*zeta*(eta + 1)*(xi - 1))/2, (zeta*(eta^2 - 1)*(xi^2 - 1))/2 + ((eta^2 - 1)*(xi^2 - 1)*(zeta - 1))/2, (zeta*(eta^2 - 1)*(xi^2 - 1))/2 + ((eta^2 - 1)*(xi^2 - 1)*(zeta + 1))/2, eta*zeta*(xi^2 - 1)*(eta - 1), eta*zeta*(xi^2 - 1)*(eta + 1), xi*zeta*(eta^2 - 1)*(xi - 1), xi*zeta*(eta^2 - 1)*(xi + 1), -2*zeta*(eta^2 - 1)*(xi^2 - 1)];
            dNe=[dNe1; dNe2; dNe3];
%%
x=[-1 0 1]; 
y=x;
z=x;
for i=1:3
    zeta=x(i);
    for j=1:3
        eta=y(i);
        for k=1:3
            xi=x(i);
            Ne=[ (eta*xi*zeta*(eta - 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(eta - 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(eta + 1)*(xi + 1)*(zeta - 1))/8, (eta*xi*zeta*(eta + 1)*(xi - 1)*(zeta - 1))/8, (eta*xi*zeta*(eta - 1)*(xi - 1)*(zeta + 1))/8, (eta*xi*zeta*(eta - 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(eta + 1)*(xi + 1)*(zeta + 1))/8, (eta*xi*zeta*(eta + 1)*(xi - 1)*(zeta + 1))/8, -(eta*zeta*(xi^2 - 1)*(eta - 1)*(zeta - 1))/4, -(xi*zeta*(eta^2 - 1)*(xi + 1)*(zeta - 1))/4, -(eta*zeta*(xi^2 - 1)*(eta + 1)*(zeta - 1))/4, -(xi*zeta*(eta^2 - 1)*(xi - 1)*(zeta - 1))/4, -(eta*zeta*(xi^2 - 1)*(eta - 1)*(zeta + 1))/4, -(xi*zeta*(eta^2 - 1)*(xi + 1)*(zeta + 1))/4, -(eta*zeta*(xi^2 - 1)*(eta + 1)*(zeta + 1))/4, -(xi*zeta*(eta^2 - 1)*(xi - 1)*(zeta + 1))/4, -(eta*xi*(zeta^2 - 1)*(eta - 1)*(xi - 1))/4, -(eta*xi*(zeta^2 - 1)*(eta - 1)*(xi + 1))/4, -(eta*xi*(zeta^2 - 1)*(eta + 1)*(xi + 1))/4, -(eta*xi*(zeta^2 - 1)*(eta + 1)*(xi - 1))/4, (zeta*(eta^2 - 1)*(xi^2 - 1)*(zeta - 1))/2, (zeta*(eta^2 - 1)*(xi^2 - 1)*(zeta + 1))/2, (eta*(xi^2 - 1)*(zeta^2 - 1)*(eta - 1))/2, (eta*(xi^2 - 1)*(zeta^2 - 1)*(eta + 1))/2, (xi*(eta^2 - 1)*(zeta^2 - 1)*(xi - 1))/2, (xi*(eta^2 - 1)*(zeta^2 - 1)*(xi + 1))/2, -(eta^2 - 1)*(xi^2 - 1)*(zeta^2 - 1)];
            sum(Ne)
        end
    end
end
%%
Nx=1;
Ny=1;
Nz=1;
IEN=quadIEN(Nx,Ny,Nz);
%%
syms e1 e2 r s Ne dNe;
Ne=[1/4*e1*e2*[(e1-1)*(e2-1) (e1+1)*(e2-1) (e1+1)*(e2+1) (e1-1)*(e2+1)] ...
            1/2*[e2*(1-e1^2)*(e2-1) e1*(e1+1)*(1-e2^2) e2*(1-e1^2)*(e2+1) ...
            e1*(e1-1)*(1-e2^2)] (1-e1^2)*(1-e2^2)];
        dNe=[1/4*e2*[(2*e1-1)*(e2-1) (2*e1+1)*(e2-1) (2*e1+1)*(e2+1) (2*e1-1)*(e2+1)] ...
            1/2*[e2*(-2*e1)*(e2-1) (2*e1+1)*(1-e2^2) e2*(-2*e1)*(e2+1) ...
            (2*e1-1)*(1-e2^2)] (-2*e1)*(1-e2^2);1/4*e1*[(e1-1)*(2*e2-1) (e1+1)*(2*e2-1) (e1+1)*(2*e2+1) (e1-1)*(2*e2+1)] ...
            1/2*[(1-e1^2)*(2*e2-1) e1*(e1+1)*(-2*e2) (1-e1^2)*(2*e2+1) ...
            e1*(e1-1)*(-2*e2)] (1-e1^2)*(-2*e2)];
        Ne=subs(Ne,e1,s);
        Ne=subs(Ne,e2,r);
        dNe=subs(dNe,e1,s);
        dNe=subs(dNe,e2,r);
%%
syms r s Ne dNe q;
Ne=[ r*(2*r - 1), s*(2*s - 1), (2*r + 2*s - 1)*(r + s - 1), 4*r*s, -4*s*(r + s - 1), -4*r*(r + s - 1)];
dNe=[ 4*r - 1, 0, 4*r + 4*s - 3, 4*s, -4*s, 4 - 4*s - 8*r; 0, 4*s - 1, 4*r + 4*s - 3, 4*r, 4 - 8*s - 4*r, -4*r];
Ne=subs(Ne,r,q);
Ne=subs(Ne,s,r);
Ne=subs(Ne,q,s);
dNe=subs(dNe,r,q);
dNe=subs(dNe,s,r);
dNe=subs(dNe,q,s)