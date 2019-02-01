%%
dt=0.5*10^-4;
Q=7*5.21;
T=3;
sat=2*T/3;
fin=4*T/3;
x=0:dt:sat;
slope=Q/sat;
y=-slope*x;
z=0:dt:fin;
B=-Q*ones(1,length(z)-length(x));
W=[y'; B']';
A=[z' W'];
fileID = fopen('cfdflow.txt','w');
fprintf(fileID,'%f %f\n',A');
fclose(fileID);