%% set up linear solve
clear all
global ID IEN nNodes nDoF EBC g Params Coord NBC; 
tol=10^-8;
elementtype='hex'; %set hex or tet - note tet performs poorly for first order
order=2; %set to 1 (linear) or 2 (quadratic) elements
rcm=0; % 1 to permute with reverse Cuthill-McKee, 0 otherwise
method='cgapprox';
offdiags=0; % set number of offdiagonals, 0 is a Jacobi preconditioner
ProblemDefinition(elementtype,order);
[K,F]=Assembly(elementtype,order);
blk=6; % sets size of blocks, 1 corresponds to Jacobi preconditioner
%K=(K+K')/2;
% permute reverse Cuthill-McKee
if rcm==1
    [R, R2]=revCM(K);
    F=F(R);
    K=K(R,R);
end
%% solve linear system
d2=K\F; %to compare
nmax=length(F);
if (strcmpi(method,'minres'))
    M=spalloc(nmax,nmax,3*nmax);
    for j=1:nmax
        M(j,:)=l_sparse_inverse(K,j,tol,nmax,offdiags);
    end
    M=(M+M')/2;
    [ d1, istop, itn, rnorm, Arnorm, Anorm, Acond, ynorm ] =  minres( K, F, M, 0, 0, 0, nmax, tol );
end
if (strcmpi(method,'cgnopre'))
    d1=zeros(nmax,1);
    [d1,niter]=conj_g(K,d1,F,nmax,tol);
end
if (strcmpi(method,'cgdiag'))
    d1=zeros(nmax,1);
    [d1,niter]=conj_diag_pre(K,d1,F,nmax,tol);
end
if (strcmpi(method,'cgapprox'))
    d1=zeros(nmax,1);
    [d1,niter]=conj_multidiag_pre(K,d1,F,nmax,tol);
end
if (strcmpi(method,'cgblock'))
    d1=zeros(nmax,1);
    [d1,niter]=conj_block(K,d1,F,nmax,tol,blk);
end
%permute back
if rcm==1
    d2=d2(R2);
    d1=d1(R2);   
end
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
c  = Params.c;
t=   Params.t;
Ix= 4/3*c^3*t;
E = Params.E;
grav=Params.grav;
andef=grav*x.^2/24/Ix/L/E.*(2*L^2+(2*L-x).^2);
figure
plot(x,defl);
hold on
plot(x, andef);
%% convergence
du=zeros(order*Nx+1,1);
switch order
    case 1
        N1e = -1/2;
        N2e =  1/2;
        for i=1:order*Nx
            J=(x(i+1)-x(i))/2;
            du(i+1)=(N2e*defl(i+1)+N1e*defl(i))/J;
        end
    case 2 % 
        N1e = -1/2;
        N2e =  1/2;
        M1e=1/2;
        M2e=-2;
        M3e=3/2;
        for i=1:order*Nx
            J=x(i+1)-x(i);
            if(mod(i,2)==1) %middle node
                du(i+1)=(N2e*defl(i+2)+N1e*defl(i))/J;
            else
                du(i+1)=(M3e*defl(i+1)+M2e*defl(i)+M1e*defl(i-1))/J;
            end
        end
    otherwise
        error('choose order 1 or 2');
end
andder=(grav*x.*((2*L - x).^2 + 2*L^2))/(12*E*Ix*L) - (grav*x.^2.*(4*L - 2*x))/(24*E*Ix*L);
figure
plot(x,du);
hold on
plot(x,andder);
dif=defl-andef;
derdif=du-andder;
L2=sqrt(dif'*dif*L/(order*Nx+1));
H1=sqrt((dif'*dif+derdif'*derdif)*L/(order*Nx+1));