function f=SampleBodyForce(X)
global C;
E=C(1,1);
v=C(1,2);
x=X(1);
y=X(2);
z=X(3);
global Params;
L = Params.L;
c = Params.c;
t = Params.t;
grav=Params.grav;
%f=grav/L/c/t/4; %change to whatever function you want
%f=[ (E*pi^2*cos((pi*x)/(2*L))*(v - 1))/(4*L^2*(2*v - 1)*(v + 1)), (E*pi^2*cos((pi*y)/(2*c))*(v - 1))/(4*c^2*(2*v - 1)*(v + 1)), (E*pi^2*cos((pi*z)/(2*t))*(v - 1))/(4*t^2*(2*v - 1)*(v + 1))];
f =[ (E*pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*L^2) - (E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*t^2) - (pi^2*sin((pi*x)/(2*L))*cos((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(4*L*t)))/2 - (E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*c^2) - (pi^2*sin((pi*x)/(2*L))*sin((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*L*c)))/2, (E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*L^2) - (pi^2*sin((pi*x)/(2*L))*sin((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*L*c)))/2 - (E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*t^2) - (pi^2*cos((pi*x)/(2*L))*sin((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(4*c*t)))/2 - (E*pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*c^2), (E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*L^2) - (pi^2*sin((pi*x)/(2*L))*cos((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(4*L*t)))/2 - (E*((pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*c^2) - (pi^2*cos((pi*x)/(2*L))*sin((pi*y)/(2*c))*sin((pi*z)/(2*t)))/(4*c*t)))/2 - (E*pi^2*cos((pi*x)/(2*L))*cos((pi*y)/(2*c))*cos((pi*z)/(2*t)))/(4*t^2)];
end