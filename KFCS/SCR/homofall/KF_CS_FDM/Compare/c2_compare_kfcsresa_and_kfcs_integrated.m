clc
clear
close all
printfigure = 0;

load f_random_simpling_KF_fdm
load f_e_scr_kf_fdm_integrated_1
load Messwerte

Dt = 0.1;

f = f(:,1:Dt/dt:end);
f1 = f_random_simpling_KF_fdm;
f2 = f_e_mt_kf_fdm_integrated_1;

NAME = 'comapre random_KF and KFCS_integrated';
l1 = 'random sampling KF';
l2 = 'KFCS integrated';




x = 0 : dx : 10;



error_1 = f - f1;
error_2 = f - f2;


error1 = zeros(size(f,2),1);
error2 = zeros(size(f,2),1);

for i = 1 : size(f,2)
    error1(i) = norm(error_1(:,i)) / norm(f(:,i));
    error2(i) = norm(error_2(:,i)) / norm(f(:,i));
end


mean_error_1 = mean(error1);
mean_error_2 = mean(error2);

t = 0 : Dt : 20;
figure
plot(t,error1,'r-','LineWidth',5)
hold on
plot(t,error2,'b-','LineWidth',5)
hold on
plot(t,ones(length(t))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_2,'b-.','LineWidth',2)
hold on
legend(l1,l2)
setplt('Error','$t$','Error',[NAME,' error'],printfigure)