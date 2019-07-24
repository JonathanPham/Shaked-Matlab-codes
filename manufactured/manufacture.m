%%
clear all;
syms x y z u u2 d_u d_u2 epsil sigma Cmat c L t E v;
u=[cos(pi*x/2/L);cos(pi*y/2/c);cos(pi*z/2/t)];
d_u=[diff(u,x),diff(u,y),diff(u,z)];
epsil=diag(d_u);
sigma1=E/(1+v)/(1-2*v)*(epsil(1)*(1-v)-v*epsil(2)-v*epsil(3));
sigma2=E/(1+v)/(1-2*v)*(epsil(2)*(1-v)-v*epsil(1)-v*epsil(3));
sigma3=E/(1+v)/(1-2*v)*(epsil(3)*(1-v)-v*epsil(2)-v*epsil(1));
f=-[diff(sigma1,x),diff(sigma2,y),diff(sigma3,z)];
h=[sigma1,sigma2,sigma3];