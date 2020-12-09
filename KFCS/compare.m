clc
clear
close all
printfigure =1;

load Messwerte

load f_e_kalman_fdm_1D_12
load f_e_kalman_fdm_1D_24
load f_e_kalman_fdm_1D_36
load f_random_sampling_KF_fdm

Dt = 0.1;

f = f(:,1:Dt/dt:end);
f1 = f_e_kalman_fdm_1D_12;
f2 = f_e_kalman_fdm_1D_24;
f3 = f_e_kalman_fdm_1D_36;
f4 = f_random_sampling_KF_fdm;


NAME = 'Compare_KF_rsKF';
l1 = 'KF FDM 12';
l2 = 'KF FDM 24';
l3 = 'KF FDM 36';
l4 = 'rsKF FDM 12';




x = 0 : dx : 10;



error_1 = f - f1;
error_2 = f - f2;
error_3 = f - f3;
error_4 = f - f4;


error1 = zeros(size(f,2),1);
error2 = zeros(size(f,2),1);
error3 = zeros(size(f,2),1);
error4 = zeros(size(f,2),1);

for i = 1 : size(f,2)
    error1(i) = norm(error_1(:,i)) / norm(f(:,i));
    error2(i) = norm(error_2(:,i)) / norm(f(:,i));
    error3(i) = norm(error_3(:,i)) / norm(f(:,i));
    error4(i) = norm(error_4(:,i)) / norm(f(:,i));
end


mean_error_1 = mean(error1);
mean_error_2 = mean(error2);
mean_error_3 = mean(error3);
mean_error_4 = mean(error4);

t = 0 : Dt : 20;
figure
plot(t,error1,'r-','LineWidth',5)
hold on
plot(t,error2,'g-','LineWidth',5)
hold on
plot(t,error3,'c-','LineWidth',5)
hold on
plot(t,error4,'k-','LineWidth',5)
hold on
plot(t,ones(length(t))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_2,'g-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_3,'c-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_4,'k-.','LineWidth',2)
legend(l1,l2,l3,l4)
setplt('','$t$','$\epsilon$',[NAME,' error'],printfigure)

close all

