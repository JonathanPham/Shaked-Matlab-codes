function f=SampleBodyForce(X)
x=X(1);
y=X(2);
z=X(3);
global Params;
L = Params.L;
c = Params.c;
t = Params.t;
grav=Params.grav;
%force=grav/L/c/t/4; %change to whatever function you want
f=[ -(pi^2*cos((pi*x)/L))/(E*L^2), -(pi^2*cos((pi*y)/c))/(E*c^2), -(pi^2*cos((pi*z)/t))/(E*t^2)];
end