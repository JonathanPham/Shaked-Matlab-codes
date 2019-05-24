%% set up linear solve
global ID IEN nNodes nDoF EBC g Params Coord grav; 
tol=1e-8;
elementtype='tet'; %set hex or tet - note tet performs poorly for first order
order=2; %set to 1 (linear) or 2 (quadratic) elements
ProblemDefinition(elementtype,order);
[K,F]=Assembly(elementtype,order);
%global IEN;
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
mmx=max(d1-d2);
%% constructing solution vector
u=zeros(nNodes,nDoF);
I=(EBC==0);
u(I)=d1(ID(I));
u(~I)=g(~I);
% post processing and comparing to analytical solution
indy=Coord(:,2)==0;
indz=Coord(:,3)==0;
ind=indy&indz;
x = Coord(ind,1);
defl=u(ind,2);
[x, ind]=sort(x);
defl=defl(ind);
Nx = Params.Nx;
L  = Params.L;
plot(x,defl);
c  = Params.c;
t=   Params.t;
Ix= 4/3*c^3*t;
E = Params.E;
andef=grav*x.^2/24/Ix/L/E.*(2*L^2+(2*L-x).^2);
hold on
%figure
plot(x, andef);
% %% convergence
% du=zeros(order*Nx+1,1);
% switch order
%     case 1
%         N1e = -1/2;
%         N2e =  1/2;
%         for i=1:order*Nx
%             J=(x(i+1)-x(i))/2;
%             du(i+1)=(N2e*defl(i+1)+N1e*defl(i))/J;
%         end
%     case 2 % this seems still accurate, but we could change it if we cared
%         N1e = -1/2;
%         N2e =  1/2;
%         for i=1:order*Nx
%             J=(x(i+1)-x(i))/2;
%             du(i+1)=(N2e*defl(i+1)+N1e*defl(i))/J;
%         end
%     otherwise
%         error('choose order 1 or 2');
% end
% andder=(grav*x.*((2*L - x).^2 + 2*L^2))/(12*E*Ix*L) - (grav*x.^2.*(4*L - 2*x))/(24*E*Ix*L);
% figure
% plot(x,du);
% hold on
% plot(x,andder);
% dif=defl'-andef;
% derdif=du'-andder;
% L2=sqrt(dif*dif'*L/(Nx+1));
% H1=sqrt((dif*dif'+derdif*derdif')*L/(Nx+1));