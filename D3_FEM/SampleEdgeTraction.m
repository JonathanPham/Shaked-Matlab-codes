function [h]=SampleEdgeTraction(x,y,z,e,k) % I think this is correct
global Params  NBC;
P = Params.P;
L = Params.L;
c = Params.c;
t = Params.t;
K1 =  3*P*L*t/(2*c^2);
K2 = -3*P*t/(4*c);
K3=0;
switch NBC(k,e)
    case 1
        h = [K1*(y/c),  K2*(1-(y/c)^2) , K3]';
        
    case 2
        h = [       0, -K2*(1-(y/c)^2), K3 ]';
        
end
end