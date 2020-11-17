clc
clear
close all
printfigure = 1;

load Messwerte_rh
load f_rh_kalman_fdm_1D_12
load f_rh_kalman_fdm_1D_24
load f_rh_random_simpling_KF_fdm
load f_rh_scr_kf_fdm_integrated_2
load f_rh_scr_kf_fdm_integrated_iterativ_2

Dt = 0.1;

f = f_sr(:,1:Dt/dt:end);
f1 = f_rh_kalman_fdm_1D_12;
f3 = f_rh_kalman_fdm_1D_24;
f2 = f_rh_random_simpling_KF_fdm;
f4 = f_rh_scr_kf_fdm_integrated_2;
f5 = f_rh_scr_kf_fdm_integrated_iterativ_2;


NAME = 'comapre noisy signal';
l1 = 'KF 12';
l3 = 'KF 24';
l2 = 'rsKF 12';
l4 = 'KFCS 12';
l5 = 'KFCS iterativ 12';




x = 0 : dx : 10;


error_1 = f - f1;
error_2 = f - f2;
error_3 = f - f3;
error_4 = f - f4;
error_5 = f - f5;


error1 = zeros(size(f,2),1);
error2 = zeros(size(f,2),1);
error3 = zeros(size(f,2),1);
error4 = zeros(size(f,2),1);
error5 = zeros(size(f,2),1);

for i = 1 : size(f,2)
    error1(i) = norm(error_1(:,i)) / norm(f(:,i));
    error2(i) = norm(error_2(:,i)) / norm(f(:,i));
    error3(i) = norm(error_3(:,i)) / norm(f(:,i));
    error4(i) = norm(error_4(:,i)) / norm(f(:,i));
    error5(i) = norm(error_5(:,i)) / norm(f(:,i));
end


mean_error_1 = mean(error1);
mean_error_2 = mean(error2);
mean_error_3 = mean(error3);
mean_error_4 = mean(error4);
mean_error_5 = mean(error5);

t = 0 : Dt : 20;
figure
plot(t,error1,'r-','LineWidth',5)
hold on
plot(t,error2,'b-','LineWidth',5)
hold on
% plot(t,error3,'g-','LineWidth',5)
% hold on
plot(t,error4,'c-','LineWidth',5)
hold on
plot(t,error5,'m-','LineWidth',5)
hold on
plot(t,ones(length(t))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_2,'b-.','LineWidth',2)
hold on
% plot(t,ones(length(t))*mean_error_3,'g-.','LineWidth',2)
% hold on
plot(t,ones(length(t))*mean_error_4,'c-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_5,'m-.','LineWidth',2)
hold on
legend(l1,l2,l4,l5)
setplt('Fehler','$t$','Fehler',[NAME,' error'],printfigure)