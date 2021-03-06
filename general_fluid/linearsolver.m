%% residual dt=10^-5
%data
clear all;
gmres8=[1.19E+03
3.17E-03
2.46E-09];
gmres6= [1.19E+03
1.36E-01
2.92E-06];
gmres4=[1.19E+03
6.01E+00
3.40E-03
5.94E-06];
gmres2=[1.19E+03
4.38E+02
1.28E+01
4.16E-01
1.60E-02
4.88E-04
3.58E-05
2.01E-06];
SIMPLE=[1.19E+03
1.06E+03
9.01E+02
7.45E+02
6.07E+02
4.88E+02
3.96E+02
3.16E+02
2.43E+02
1.92E+02
1.52E+02
1.15E+02
9.21E+01
7.06E+01
5.62E+01
4.25E+01
3.35E+01
2.67E+01
2.02E+01
1.66E+01
1.23E+01];
asm=[1.19E+03
9.03E+02
6.64E+02
6.16E+02
5.60E+02
5.19E+02
4.86E+02
4.54E+02
4.33E+02
4.15E+02
4.00E+02
3.85E+02
3.71E+02
3.59E+02
3.49E+02
3.40E+02
3.31E+02
3.23E+02
3.16E+02
3.09E+02
3.02E+02];
times=[50.5	29.05	20.95	13.55	76.5	51];

%% residual dt=10^-1
gmres8=[1.59E+03
2.32E-02
3.11E-10];
gmres6=[1.59E+03
1.41E+00
1.29E-06];
gmres4=[1.59E+03
1.37E+02
1.51E-02
1.96E-06];
gmres2=[1.59E+03
4.42E+02
1.32E+02
2.34E+01
7.30E-01
8.77E-02
1.07E-02
1.64E-03
1.98E-04
1.01E-05];
SIMPLE=[1.59E+03
1.53E+03
8.93E+02
6.31E+02
3.81E+02
2.06E+02
5.78E+01
2.25E+01
8.69E+00
4.21E+00
1.78E+00
1.05E+00
2.90E-01
7.95E-02
3.27E-02
1.46E-02
3.93E-03
1.09E-03
3.46E-04
6.81E-05
2.30E-05];
asm=[1.59E+03
9.42E+02
7.85E+02
6.80E+02
6.18E+02
5.66E+02
5.29E+02
4.97E+02
4.70E+02
4.49E+02
4.29E+02
4.12E+02
3.98E+02
3.85E+02
3.73E+02
3.62E+02
3.53E+02
3.44E+02
3.36E+02
3.29E+02
3.15E+02];
times=[140.5	100.5	89	81	174	64];
%% analysis

%normalize
asm=asm/asm(1);
SIMPLE=SIMPLE/SIMPLE(1);
gmres2=gmres2/gmres2(1);
gmres4=gmres4/gmres4(1);
gmres6=gmres6/gmres6(1);
gmres8=gmres8/gmres8(1);
%plot
thresh=10^-8*ones(1,length(asm));
figure;
hold on
set(gca, 'YScale', 'log')
plot(0:length(gmres8)-1,gmres8,'-*');
plot(0:length(gmres6)-1,gmres6,'-o');
plot(0:length(gmres4)-1,gmres4,'-s');
plot(0:length(gmres2)-1,gmres2,'-d');
plot(0:length(SIMPLE)-1,SIMPLE,'-X');
plot(0:length(asm)-1,asm,'-p');
plot(0:length(asm)-1,thresh,'--');
legend(strcat('\delta=10^{-8}, t=',num2str(times(1))),strcat('\delta=10^{-6}, t=',num2str(times(2))),strcat('\delta=10^{-4}, t=',num2str(times(3))),strcat('\delta=10^{-2}, t=',num2str(times(4))),strcat('SIMPLE \delta=10^{-8}, t=',num2str(times(5))),strcat('Additive Schwarz, t=',num2str(times(6))),'convergence','Location','southeast');
title('Convergence for dt=10^{-5}');
xlabel('iteration');
ylabel('||r||_2/||r_0||_2');