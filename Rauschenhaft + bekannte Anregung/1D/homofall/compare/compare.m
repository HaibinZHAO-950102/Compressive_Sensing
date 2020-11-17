clc
clear
close all
printfigure = 1;

load Messwerte_rh
load f_rh_kalman_fdm_1D_12
load f_rh_scr_kf_fdm_integrated_iterativ_2

Dt = 1;

f = f_sr(:,1:Dt/dt:end);
f1 = f_rh_kalman_fdm_1D_12;
f5 = f_rh_scr_kf_fdm_integrated_iterativ_2;


NAME = 'comapre noisy signal';
l1 = 'KF 12';
l5 = 'KFCS 12';




x = 0 : dx : 10;


error_1 = f - f1;
error_5 = f - f5;


error1 = zeros(size(f,2),1);
error5 = zeros(size(f,2),1);

for i = 1 : size(f,2)
    error1(i) = norm(error_1(:,i)) / norm(f(:,i));
    error5(i) = norm(error_5(:,i)) / norm(f(:,i));
end


mean_error_1 = mean(error1);
mean_error_5 = mean(error5);

t = 0 : Dt : 200;
figure
plot(t,error1,'r-','LineWidth',5)
hold on
plot(t,error5,'b-','LineWidth',5)
hold on
plot(t,ones(length(t))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_5,'b-.','LineWidth',2)
hold on
legend(l1,l5)
setplt('Fehler','$t$','Fehler',[NAME,' error'],printfigure)