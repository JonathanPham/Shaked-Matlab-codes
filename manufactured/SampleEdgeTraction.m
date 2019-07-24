function [h]=SampleEdgeTraction(x,y,z,e,k) % I think this is correct
global Params  NBC;
P = Params.P;
L = Params.L;
c = Params.c;
t = Params.t;
K1 =  P/L/c/2;
K2 = -P/L/c/2;
switch NBC(e,k)
    case 1
        h = -pi/2*[       0, sin(pi*y/2/c)/c, sin(pi*z/2/t)/t]';
    case 3
        h = -pi/2*[sin(pi*x/2/L)/L, 0 , sin(pi*z/2/t)/t]';
        
    case 5
        h = -pi/2*[sin(pi*x/2/L)/L, sin(pi*y/2/c)/c, 0]';
        
end
end