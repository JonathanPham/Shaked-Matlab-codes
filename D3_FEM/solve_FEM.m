%% set up linear solve
global IEN ID nNodes nDoF EBC g Params Coord grav; 
tol=1e-12;
ProblemDefinition();
[K,F]=Assembly();
%% solve linear system
d2=K\F; %to compare
nmax=length(F);
% M=spalloc(nmax,nmax,3*nmax);
% for j=1:nmax
%     M(j,:)=l_sparse_inverse(K,j,tol,nmax,1);
% end
% guess=M*F;
% [d1]=gmres(M*K,guess,guess,nmax,tol);
d1=zeros(nmax,1);
d1=conj_g(K,d1,F,nmax,tol);

max(d1-d2);
%% constructing solution vector
u=zeros(nNodes,nDoF);
I=(EBC==0);
u(I)=d1(ID(I));
u(~I)=g(~I);
% post processing and comparing to analytical solution
indy=Coord(:,2)==0;
indz=Coord(:,3)==0;
ind=indy.*indz;
defl=u(ind==1,2);
Nx = Params.Nx;
L  = Params.L;
x = linspace(0,L,Nx+1);
plot(x,defl);
c  = Params.c;
t=   Params.t;
Ix= 4/3*c^3*t;
E = Params.E;
andef=grav*x.^2/24/Ix/L/E.*(2*L^2+(2*L-x).^2);
hold on
%figure
plot(x, andef);
%% convergence
N1e = -1/2;
N2e =  1/2;
du=zeros(Nx+1,1);
for i=1:Nx
    J=(x(i+1)-x(i))/2;
    du(i+1)=(N2e*defl(i+1)+N1e*defl(i))/J;
end
andder=(grav*x.*((2*L - x).^2 + 2*L^2))/(12*E*Ix*L) - (grav*x.^2.*(4*L - 2*x))/(24*E*Ix*L);
figure 
plot(x,du);
hold on
plot(x,andder);
%% norm calculation
dif=defl'-andef;
derdif=du'-andder;
L2=sqrt(dif*dif'*L/(Nx+1));
H1=sqrt((dif*dif'+derdif*derdif')*L/(Nx+1));