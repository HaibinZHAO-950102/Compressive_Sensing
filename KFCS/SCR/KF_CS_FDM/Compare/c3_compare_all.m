clc
clear
close all
printfigure = 1;

load f_random_simpling_KF_fdm
load f_e_scr_kf_fdm
load f_e_scr_kf_fdm_integrated_1
load f_e_kalman_fdm_1D_12
load f_e_kalman_fdm_1D_40
load Messwerte

Dt = 0.1;

f = f(:,1:Dt/dt:end);
f2 = f_random_simpling_KF_fdm;
f5 = f_e_mt_kf_fdm_integrated_1;
f1 = f_e_kalman_fdm_1D_12;
f4 = f_e_kalman_fdm_1D_40;
f3 = f_e_scr_kf_fdm;


NAME = 'comapre all';
l2 = 'random sampling KF';
l5 = 'KFCS integrated';
l1 = 'KF 12';
l4 = 'KF 40';
l3 = 'KFCS';



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
plot(t,error3,'g-','LineWidth',5)
hold on
plot(t,error4,'c-','LineWidth',5)
hold on
plot(t,error5,'m-','LineWidth',5)
hold on
plot(t,ones(length(t))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_2,'b-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_3,'g-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_4,'c-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_5,'m-.','LineWidth',2)
hold on
legend(l1,l2,l3,l4,l5)
setplt('Error','$t$','Error',[NAME,' error'],printfigure)
