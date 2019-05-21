%% parameters
clear all;
N=1000;
d=zeros(1,N+1);
d(1)=1;
v=zeros(1,N+1);
v(1)=1;
k=pi^2;
a=zeros(1,N+1);
a(1)=-k*d(1);
roinf=0.8;
tn=0.4;
dt=tn/N;
t=0:dt:tn;
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
%% all
gamma=1/2-am+af;
beta=1/4*(1-am+af)^2;
for i=1:N
    d(i+1)=((1+dt^2*beta*af*k/(am-1))*d(i)+dt*v(i)+(0.5-beta+beta*am/(am-1))*dt^2*a(i))/(1+(1-af)/(1-am)*k*beta*dt^2);
    a(i+1)=(am*a(i)+k*((1-af)*d(i+1)+af*d(i)))/(am-1);
    v(i+1)=v(i)+dt*((1-gamma)*a(i)+gamma*a(i+1));
end
%% analytic
dan=1/pi*sin(pi*t)+cos(pi*t);
logdiff=log10(abs(d-dan));
plot(t(101:1001),logdiff(101:1001));