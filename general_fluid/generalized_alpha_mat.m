%% parameters
clear all;
N=1000;
tn=0.4;
dt=tn/N;
t=0:dt:tn;
d=zeros(1,N+1);
d(1)=1;
vdt=zeros(1,N+1);
vdt(1)=1*dt;
k=pi^2;
adt2=zeros(1,N+1);
adt2(1)=-k*d(1)*dt^2;
roinf=0.8;
figure 
hold on
%% HHT
am=0;
af=(1-roinf)/(1+roinf);
%% WBZ
af=0;
am=(roinf-1)/(1+roinf);
%% gen alpha
am=(2*roinf-1)/(1+roinf);
af=(roinf)/(1+roinf);
%% legacy
d(i+1)=((1+dt^2*beta*af*k/(am-1))*d(i)+dt*v(i)+(0.5-beta+beta*am/(am-1))*dt^2*a(i))/(1+(1-af)/(1-am)*k*beta*dt^2);
a(i+1)=(am*a(i)+k*((1-af)*d(i+1)+af*d(i)))/(am-1);
v(i+1)=v(i)+dt*((1-gamma)*a(i)+gamma*a(i+1));
%% matrix method
gamma=1/2-am+af;
beta=1/4*(1-am+af)^2;
mat=zeros(3);
facx=1+(1-af)/(1-am)*k*beta*dt^2;
mat(1,1)=(1+dt^2*beta*af*k/(am-1))/facx;
mat(1,2)=1/facx;
mat(1,3)=(0.5+beta/(am-1))/facx;
mat(3,1)=dt^2*k/(am-1)*(af+(1-af)*mat(1,1));
mat(3,2)=dt^2*(1-af)*k/(am-1)*mat(1,2);
mat(3,3)=dt^2*(1-af)*k/(am-1)*mat(1,3)+am/(1-am);
mat(2,1)=gamma*mat(3,1);
mat(2,2)=1+gamma*mat(3,2);
mat(2,3)=1-gamma+gamma*mat(3,3);
for i=1:N
    vec=mat*[d(i);vdt(i);adt2(i)];
    d(i+1)=vec(1);
    vdt(i+1)=vec(2);
    adt2(i+1)=vec(3);   
end
%% analytic
dan=1/pi*sin(pi*t)+cos(pi*t);
logdiff=log10(abs(d-dan));
plot(t(101:1001),logdiff(101:1001));
%% spectrum
clear all;
N=1000;
tn=2;
dt=tn/N;
roinf=0.8;
k=pi^2;
T=2*pi/sqrt(k);
spect=zeros(1,N);
am=(2*roinf-1)/(1+roinf);
af=(roinf)/(1+roinf);
gamma=1/2-am+af;
beta=1/4*(1-am+af)^2;
mat=zeros(3);
for i=1:N
    k=pi^2*delt(i);
    facx=1+(1-af)/(1-am)*k*beta*dt^2;
    mat(1,1)=(1+dt^2*beta*af*k/(am-1))/facx;
    mat(1,2)=1/facx;
    mat(1,3)=(0.5+beta/(am-1))/facx;
    mat(3,1)=dt^2*k/(am-1)*(af+(1-af)*mat(1,1));
    mat(3,2)=dt^2*(1-af)*k/(am-1)*mat(1,2);
    mat(3,3)=dt^2*(1-af)*k/(am-1)*mat(1,3)+am/(1-am);
    mat(2,1)=gamma*mat(3,1);
    mat(2,2)=1+gamma*mat(3,2);
    mat(2,3)=1-gamma+gamma*mat(3,3);
    spect(i)=max(abs(eig(mat)));
end