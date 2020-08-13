clc
clear
close all
printfigure = 0;

load Messwerte_rh
load f_rh_kalman_fdm_1D_12
load f_rh_kalman_fdm_1D_24
load f_rh_random_simpling_KF_fdm
load f_rh_scr_kf_fdm_integrated_2
load f_rh_scr_kf_fdm_integrated_iterativ_2

Dt = 0.1;

f = f_sr(:,1:Dt/dt:end);
f1 = f_rh_kalman_fdm_1D_12;
f2 = f_rh_kalman_fdm_1D_24;
f3 = f_rh_random_simpling_KF_fdm;
f4 = f_rh_scr_kf_fdm_integrated_2;
f5 = f_rh_scr_kf_fdm_integrated_iterativ_2;


NAME = 'comapre noisy signal';
l1 = 'KF 12';
l2 = 'KF 24';
l3 = 'ramdom sampling KF';
l4 = 'KFCS';
l5 = 'KFCS iterativ';




x = 0 : dx : 10;


Error_1 = f - f1;
Error_2 = f - f2;
Error_3 = f - f3;
Error_4 = f - f4;
Error_5 = f - f5;



for i = 1 : size(f,2)
    error_1(i) = norm(Error_1(:,i));
    error_2(i) = norm(Error_2(:,i));
    error_3(i) = norm(Error_3(:,i));
    error_4(i) = norm(Error_4(:,i));
    error_5(i) = norm(Error_5(:,i));
end

% error_1 = 2 * (1-similariy_1);
% error_2 = 2 * (1-similariy_2);
% error_3 = 2 * (1-similariy_3);
% error_4 = 2 * (1-similariy_4);
% error_5 = 2 * (1-similariy_5);
% 
mean_error_1 = mean(error_1(25:end));
mean_error_2 = mean(error_2(25:end));
mean_error_3 = mean(error_3(25:end));
mean_error_4 = mean(error_4(25:end));
mean_error_5 = mean(error_5(25:end));

t = 0 : Dt : 20;
figure
plot(t(25:end),error_1(25:end),'r-','LineWidth',5)
hold on
plot(t(25:end),error_2(25:end),'b-','LineWidth',5)
hold on
plot(t(25:end),error_3(25:end),'g-','LineWidth',5)
hold on
plot(t(25:end),error_4(25:end),'c-','LineWidth',5)
hold on
plot(t(25:end),error_5(25:end),'k-','LineWidth',5)
hold on
plot(t(25:end),ones(length(t(25:end)))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t(25:end),ones(length(t(25:end)))*mean_error_2,'b-.','LineWidth',2)
hold on
plot(t(25:end),ones(length(t(25:end)))*mean_error_3,'g-.','LineWidth',2)
hold on
plot(t(25:end),ones(length(t(25:end)))*mean_error_4,'c-.','LineWidth',2)
hold on
plot(t(25:end),ones(length(t(25:end)))*mean_error_5,'k-.','LineWidth',2)
hold on
legend(l1,l2,l3,l4,l5)
setplt('Error','$t$','Error',[NAME,' error'],printfigure)
