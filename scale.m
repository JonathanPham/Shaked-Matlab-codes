%% data
clear all;
gmresCPU=[854
2020
4455];
gmreslin=[786
1929
4303];
asmCPU=[900,4136, nan];
asmlin=[830, 4044, nan];
simpleCPU=[939
2632
5987];
simplelin=[870
2541
5867];
elem=[48,96,192];
%% plot
figure; 
hold on;
plot(elem,gmresCPU);
plot(elem,gmreslin);
plot(elem,asmCPU);
plot(elem,asmlin);
plot(elem,simpleCPU);
plot(elem,simplelin);
legend('gmres CPU','gmres linear solve', 'a.s. CPU', 'a.s. linear', 'SIMPLE CPU','SIMPLE linear','Location','southeast');
ylabel('time');
xlabel('elements/dimension');

%% % plots
elem=[48,96,192];
timeg=[21.24324324
50.76315789
107.575];
times=[23.51351351
66.86842105
146.675];
timeasm=[22.43243243, 109.2972973, nan];

ling=100*[0.9203747073
0.954950495
0.9658810325];
linas=100*[0.9222222222, 0.9777562863, nan];
lins=100*[0.9265175719
0.9654255319
0.9799565726];
figure;
hold on;
plot(elem,timeg);
plot(elem,timeasm);
plot(elem,times);
xlabel('elements/dimension');
ylabel('time/linear solve');
legend('gmres','Additive Schwarz','SIMPLE','Location','southeast');
figure;
hold on;
plot(elem,ling);
plot(elem,linas);
plot(elem,lins);
xlabel('elements/dimension');
ylabel('% linear solve');
legend('gmres','Additive Schwarz','SIMPLE','Location','southeast');
%% ref data
hold on; 
elem=[48,96,192];
timeg2=[19.7 50.2 108];
timeasm2=[20.4 106.8 305.4];
times2=[21.6 63.6 153.6];
plot(elem,timeg2);
plot(elem,timeasm2);
plot(elem,times2);