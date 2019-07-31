function f=SampleBodyForce(X)
global Params C;
L = Params.L;
c = Params.c;
t = Params.t;
E=C(1,1);
x=X(1);
y=X(2);
z=X(3);
grav=Params.grav;
% force=grav/L/c/t/4; %change to whatever function you want
% f=[0 force 0];
f=[ E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*c^2) - (pi^2*sin((pi*x)/(2*L))*sin((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*L*c)) + E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*t^2) - (pi^2*sin((pi*x)/(2*L))*cos((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(8*L*t)) + (E*pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*L^2), E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*L^2) - (pi^2*sin((pi*x)/(2*L))*sin((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*L*c)) + E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*t^2) - (pi^2*cos((pi*x)/(2*L))*sin((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(8*c*t)) + (E*pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*c^2), E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*L^2) - (pi^2*sin((pi*x)/(2*L))*cos((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(8*L*t)) + E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(8*c^2) - (pi^2*cos((pi*x)/(2*L))*sin((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(8*c*t)) + (E*pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*t^2)];
end