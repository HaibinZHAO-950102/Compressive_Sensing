clc
clear
close all
printfigure = 1;

load Messwerte_rh
load f_e_kf
load f_e_kfcs


f = f;
f0 = zeros(size(f));
f1 = f_e_kf;
f2 = f_e_kfcs;


NAME = 'comapre KF and KFCS';
l1 = '2D KF';
l2 = '2D KFCS';


x = 1 : 128^2;

error0 = zeros(size(f,2),1);
error1 = zeros(size(f,2),1);
error2 = zeros(size(f,2),1);

for i = 1 : size(f,2)
    error0(i) = norm(f(:,i) - f0(:,i));
    error1(i) = norm(f(:,i) - f1(:,i));
    error2(i) = norm(f(:,i) - f2(:,i));
end

mean_error_0 = mean(error0);


error1 = error1 / mean_error_0;
error2 = error2 / mean_error_0;

mean_error_1 = mean(error1);
mean_error_2 = mean(error2);

t = 0 : 0.1 : 20;
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
