%%
mat=load('vpavg.csv');
velc=[mat(:,9:16:end-16),mat(:,4)];
mvelc=mean(velc,2);
umean=velmult*4.61;
rho=1.056;
fac=0.5*rho*umean^2*100;
hold on
grid on
plot(mat(:,end)/100,mvelc/umean,'g','Linewidth',2);
%%
prc=[mat(:,6:16:end-20)];
mprc=mean(prc,2);
p0=mprc(501,1);
plot(mat(:,end)/100,(mprc-p0)/fac,'g','Linewidth',2);
hold on
%%
mat=load('cut4avg.csv');
velc=[mat(:,9:16:end-16),mat(:,4)];
mvelc=mean(velc,2);
umean=velmult*4.61;
rho=1.056;
fac=0.5*rho*umean^2*100;
hold on
plot(mat(:,end-1)/100,mvelc/umean,'g','Linewidth',2);
%legend('SV time average','data1','data2','data3','data4','data5','SV refined');