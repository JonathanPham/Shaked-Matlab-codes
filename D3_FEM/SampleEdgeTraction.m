function [h]=SampleEdgeTraction(x,y,z,e,k) % I think this is correct
global Params  NBC;
P = Params.P;
L = Params.L;
c = Params.c;
t = Params.t;
K1 =  P/L/c/2;
K2 = -P/L/c/2;
switch NBC(e,k)
    case 5
        h = [0, K1 , 0]';
        
    case 6
        h = [       0, K2, 0]';
        
end
end