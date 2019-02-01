clear all;
velmult=7;
%%
mat=load('vpavg.csv');
umean=velmult*4.61;
p0=mat(501,1);
rho=1.056;
fac=0.5*rho*umean^2*100;
%figure 
plot(mat(:,23)/100,mat(:,4)/umean,'k','Linewidth',2);
xlabel('z-axis');
ylabel('normalized velocity');
hold on
%%
zvel=load('z_vel.txt');
plot(zvel(:,1),zvel(:,2)*100/umean,'Linewidth',2);
%%
grid on
xlim([-0.08 0.08]);
legend('SV','data1','data2','data3','data4','data5');
%%
grid on
xlim([-0.08 0.08]);
legend('SV coarse','SV fine','data1','data2','data3','data4','data5');
%%
plot(mat(:,23)/100,(mat(:,1)-p0)/fac,'k','Linewidth',2);
xlabel('z-axis');
ylabel('normalized pressure difference');
hold on
%%
zpr=load('z_press.txt');
plot(zpr(:,1),zpr(:,2)*10/fac,'Linewidth',2);
%%
legend('SV','SV_{refined}','data1','data2','data3','data4','data5');
%%
mat1=load('vel_cut1.csv');
mat2=load('vel_cut2.csv');
mat3=load('vel_cut3.csv');
mat4=load('vel_cut4.csv');
umean=velmult*4.61;
rd=0.006;
%%
plot(mat1(:,20)/100-rd,mat1(:,4)/umean,'k','Linewidth',2);
hold on
%plot(matb1(:,20)-rd,matb1(:,4)/umean);
xlabel('radial axis');
ylabel('u_{norm} @ z=-0.064');
xlim([-rd rd]);
%%
v1=load('cut1.txt');
plot(v1(:,1),v1(:,2)*100/umean,'Linewidth',2);
%%
grid on
legend('SV','data1','data2','data3','data4','data5');
%%
plot(mat2(:,20)/100-rd,mat2(:,4)/umean,'k','Linewidth',2);
hold on
%plot(matb2(:,20)-rd,matb2(:,4)/umean);
xlabel('radial axis');
ylabel('u_{norm} @ z=-0.008');
xlim([-rd rd]);
%%
v2=load('cut2.txt');
plot(v2(:,1),v2(:,2)*100/umean,'Linewidth',2);
%%
grid on
legend('SV coarse','SV fine','data1','data2','data3','data4','data5');
%%
plot(mat3(:,20)/100-rd,mat3(:,4)/umean,'k','Linewidth',2);
hold on
%plot(matb3(:,20)-rd,matb3(:,4)/umean);
xlabel('radial axis');
ylabel('u_{norm} @ z=0.016');
xlim([-rd rd]);
%%
v3=load('cut3.txt');
plot(v3(:,1),v3(:,2)*100/umean,'Linewidth',2);
%%
grid on
legend('SV coarse','SV fine','data1','data2','data4','data5');
%%
plot(mat4(:,20)/100-rd,mat4(:,4)/umean,'k','Linewidth',2);
hold on
%plot(matb4(:,20)-rd,matb4(:,4)/umean);
xlabel('radial axis');
ylabel('u_{norm} @ z=0.06');
xlim([-rd rd]);
%%
v4=load('cut4.txt');
plot(v4(:,1),v4(:,2)*100/umean,'Linewidth',2);
%%
legend('SV','SV_{refined}','data1','data2','data3','data4','data5');
