function [h]=SampleEdgeTraction(x,y,z,e,k) % I think this is correct
global Params  NBC;
P = Params.P;
L = Params.L;
c = Params.c;
t = Params.t;
E = Params.E;
v = Params.v;
% switch NBC(e,k)
%     case 1
%         %h = -pi/2*[       0, sin(pi*y/2/c)/c, sin(pi*z/2/t)/t]';
%         h=[ -(E*((v*pi*sin((pi*y)/(2*c)))/(2*c) + (v*pi*sin((pi*z)/(2*t)))/(2*t) + (pi*sin((pi*x)/(2*L))*(v - 1))/(2*L)))/((2*v - 1)*(v + 1)), -(E*((v*pi*sin((pi*x)/(2*L)))/(2*L) + (v*pi*sin((pi*z)/(2*t)))/(2*t) + (pi*sin((pi*y)/(2*c))*(v - 1))/(2*c)))/((2*v - 1)*(v + 1)), -(E*((v*pi*sin((pi*x)/(2*L)))/(2*L) + (v*pi*sin((pi*y)/(2*c)))/(2*c) + (pi*sin((pi*z)/(2*t))*(v - 1))/(2*t)))/((2*v - 1)*(v + 1))]';
%     case 3
%         %h = -pi/2*[sin(pi*x/2/L)/L, 0 , sin(pi*z/2/t)/t]';
%         h=[ -(E*((v*pi*sin((pi*y)/(2*c)))/(2*c) + (v*pi*sin((pi*z)/(2*t)))/(2*t) + (pi*sin((pi*x)/(2*L))*(v - 1))/(2*L)))/((2*v - 1)*(v + 1)), -(E*((v*pi*sin((pi*x)/(2*L)))/(2*L) + (v*pi*sin((pi*z)/(2*t)))/(2*t) + (pi*sin((pi*y)/(2*c))*(v - 1))/(2*c)))/((2*v - 1)*(v + 1)), -(E*((v*pi*sin((pi*x)/(2*L)))/(2*L) + (v*pi*sin((pi*y)/(2*c)))/(2*c) + (pi*sin((pi*z)/(2*t))*(v - 1))/(2*t)))/((2*v - 1)*(v + 1))]';
%     case 5
%         %h = -pi/2*[sin(pi*x/2/L)/L, sin(pi*y/2/c)/c, 0]';
%         h=[ -(E*((v*pi*sin((pi*y)/(2*c)))/(2*c) + (v*pi*sin((pi*z)/(2*t)))/(2*t) + (pi*sin((pi*x)/(2*L))*(v - 1))/(2*L)))/((2*v - 1)*(v + 1)), -(E*((v*pi*sin((pi*x)/(2*L)))/(2*L) + (v*pi*sin((pi*z)/(2*t)))/(2*t) + (pi*sin((pi*y)/(2*c))*(v - 1))/(2*c)))/((2*v - 1)*(v + 1)), -(E*((v*pi*sin((pi*x)/(2*L)))/(2*L) + (v*pi*sin((pi*y)/(2*c)))/(2*c) + (pi*sin((pi*z)/(2*t))*(v - 1))/(2*t)))/((2*v - 1)*(v + 1))]';
% end
if mod(NBC(e,k),2)==1
    h=[ -(E*pi*sin((pi*x)/(2*L)))/(2*L), -(E*pi*sin((pi*y)/(2*c)))/(2*c), -(E*pi*sin((pi*z)/(2*t)))/(2*t)]';
end
end