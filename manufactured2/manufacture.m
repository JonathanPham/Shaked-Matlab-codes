%%
clear all;
syms x y z u d_u epsil veps Cmat c L t E v;
u=[cos(pi*x/2/L),cos(pi*y/2/c),cos(pi*z/2/t)];
d_u=[diff(u(1),x),diff(u(1),y),diff(u(1),z);diff(u(2),x),diff(u(2),y),diff(u(2),z);diff(u(3),x),diff(u(3),y),diff(u(3),z)];
for i=1:3
    for j=1:3
        epsil(i,j)=1/2*(d_u(i,j)+d_u(j,i));
    end
end
veps=[epsil(1,1); epsil(2,2); epsil(3,3); 2*epsil(1,2); 2*epsil(2,3); 2*epsil(3,1)];
C(1,1)=1-v;
C(2,2)=1-v;
C(3,3)=1-v;
C(4,4)=(1-2*v)/2;
C(5,5)=(1-2*v)/2;
C(6,6)=(1-2*v)/2;
C(1,2)=v;
C(1,3)=v;
C(2,1)=v;
C(2,3)=v;
C(3,1)=v;
C(3,2)=v;
C=E/(v+1)/(1-2*v)*C;
C=subs(C,v,0);
%%
syms vsig h;
vsig=C*veps;
% sigma(1,1)=vsig(1);
% sigma(2,2)=vsig(2);
% sigma(3,3)=vsig(3);
% sigma(1,2)=vsig(4);
% sigma(2,1)=vsig(4);
% sigma(2,3)=vsig(5);
% sigma(3,2)=vsig(5);
% sigma(3,1)=vsig(6);
% sigma(1,3)=vsig(6);
% h=sum(sigma);
h(1)=vsig(1)+vsig(4)+vsig(6);
h(2)=vsig(4)+vsig(2)+vsig(5);
h(3)=vsig(6)+vsig(5)+vsig(3);
%%
syms f;
f(1)=-diff(vsig(1),x)+diff(vsig(4),y)+diff(vsig(6),z);
f(2)=-diff(vsig(4),x)+diff(vsig(2),y)+diff(vsig(5),z);
f(3)=-diff(vsig(6),x)+diff(vsig(5),y)+diff(vsig(3),z);