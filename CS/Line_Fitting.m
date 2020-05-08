clc
clear
close all

printfigure = 0;

xd_1 = sort(rand(15,1)*10);
yd_1 = 5 + xd_1 + randn(length(xd_1),1);

f1 = @(p,x,y) sum(abs(p(1)*x+p(2)-y));
fun1 = @(p)f1(p,xd_1,yd_1);

f2 = @(p,x,y) (p(1)*x+p(2)-y)'*(p(1)*x+p(2)-y);
fun2 = @(p)f2(p,xd_1,yd_1);

p1_1 = fminsearch(fun1,[1,1]);
p2_1 = fminsearch(fun2,[1,1]);

x = 0 : 0.001 : 10;
y1_1 = p1_1(1) * x + p1_1(2);
y2_1 = p2_1(1) * x + p2_1(2);

plot(xd_1,yd_1,'k.','Markersize',40)
hold on
plot(x,y1_1,'b-','LineWidth',2)
hold on
plot(x,y2_1,'r-','LineWidth',2)
legend('Data','L1 fit','L2 fit')
xlim([0 10])
ylim([0 20])
setplt('Data Fitting 1','$x$','$y$','Data_Fitting_1',printfigure)

xd_2 = [xd_1 ; 1 ; 9];
yd_2 = [yd_1 ; 18 ; 2];

fun1 = @(p)f1(p,xd_2,yd_2);
fun2 = @(p)f2(p,xd_2,yd_2);

p1_2 = fminsearch(fun1,[1,1]);
p2_2 = fminsearch(fun2,[1,1]);

y1_2 = p1_2(1) * x + p1_2(2);
y2_2 = p2_2(1) * x + p2_2(2);

figure
plot(xd_2,yd_2,'k.','Markersize',40)
hold on
plot(x,y1_2,'b-','LineWidth',2)
hold on
plot(x,y2_2,'r-','LineWidth',2)
legend('Data','L1 fit','L2 fit')
xlim([0 10])
ylim([0 20])
setplt('Data Fitting 2','$x$','$y$','Data_Fitting_2',printfigure)


