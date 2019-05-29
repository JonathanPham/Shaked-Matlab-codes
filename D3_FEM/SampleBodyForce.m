function f=SampleBodyForce(X)
global Params;
L = Params.L;
c = Params.c;
t = Params.t;
grav=Params.grav;
force=grav/L/c/t/4; %change to whatever function you want
f=[0 force 0];
end