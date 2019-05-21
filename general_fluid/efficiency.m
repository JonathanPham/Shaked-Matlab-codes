%% data
clear all
vec=1:7;
proc=(2*ones(1,7)).^vec;
CPU=[22220
11550
5781
3468
2185
1159
640.6];
lin=[19034
9934.1
4951.5
3048.6
1941.2
1018.4
567.8];
CPU2=[24900 12600 6480 3430 1870 956 449];
lin2=[21600 10900 5580 2960 1620 837 384];
CPU=CPU/CPU(1);
lin=lin/lin(1);
CPU2=CPU2/CPU2(1);
lin2=lin2/lin2(1);
ideal=2.^(-1*[0:6]);
plot(proc,CPU,'Linewidth',2);
hold on;
plot(proc,CPU2,'Linewidth',2);
plot(proc,ideal,'Linewidth',2);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log')
legend('results','reference','ideal');
xlabel('Processors');
ylabel('normalized time');
%% plot scale
figure
hold on
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log')
plot(proc,CPU);
plot(proc,lin);
plot(proc, CPU2);
plot(proc,lin2);
xlabel('Processors');
ylabel('time');
legend('total CPU', 'linear solver', 'total CPU - reference', 'linear solver - reference');
title('scaling with dt=10^{-5}');
%% plot efficiency
eff=[1
0.9619047619
0.9609064176
0.800893887
0.635583524
0.5991156169
0.5419723697];
semilogx(proc,eff);
hold on;
eff2=[1 .99 .96 .91 .83 .81 .87];
semilogx(proc,eff2);
xlabel('Processors');
ylabel('efficiency');
legend('results', 'reference');
%%
