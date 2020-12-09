clc
clear
close all
printfigure = 1;

load Messwerte_rh
load f_e_kf_25
load f_e_kf_64
load f_e_kfcs


f = f;
f0 = zeros(size(f));
f1 = f_e_kf_25;
f2 = f_e_kf_64;
f3 = f_e_kfcs;


NAME = 'comapre KF and KFCS';
l1 = '2D KF mit 25 Sensoren';
l2 = '2D KF 64';
l3 = '2D KFCS mit 25 Messungen';


x = 1 : 128^2;

error0 = zeros(size(f,2),1);
error1 = zeros(size(f,2),1);
error2 = zeros(size(f,2),1);
error3 = zeros(size(f,2),1);

for i = 1 : size(f,2)
    error0(i) = norm(f(:,i) - f0(:,i));
    error1(i) = norm(f(:,i) - f1(:,i));
    error2(i) = norm(f(:,i) - f2(:,i));
    error3(i) = norm(f(:,i) - f3(:,i));
end

mean_error_0 = mean(error0);


error1 = error1 / mean_error_0;
error2 = error2 / mean_error_0;
error3 = error3 / mean_error_0;

mean_error_1 = mean(error1);
mean_error_2 = mean(error2);
mean_error_3 = mean(error3);

t = 0 : 0.1 : 20;
figure
plot(t,error1,'r-','LineWidth',5)
hold on
plot(t,error3,'k-','LineWidth',5)
hold on
plot(t,ones(length(t))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t,ones(length(t))*mean_error_3,'k-.','LineWidth',2)
hold on
legend(l1,l3)
setplt('','$t$','$\epsilon$',[NAME,' error'],printfigure)
 close all