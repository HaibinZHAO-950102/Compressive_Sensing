clc
clear
close all

printfigure = 1;

load('Messwerte.mat')


y_1 = f(p_index,2001-50);
y_2 = f(p_index,2001);
x = 0 : dx : 10;
xm = x(p_index);
N = length(y_1);
K = 10;

WT = DWT(64,'haar');
THETA = WT^-1;

z_1 = WT * y_1;
z_2 = WT * y_2;

plot(z_1,'k.','markersize',20)
xlim([1 N])
txt = ['$t = 19.5$'];
T = text(50,8,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients of Temperature with WT $t_{1}$','$z$','$Value$','Coeffcients of Temperature with WT t1',printfigure)

figure
plot(z_2,'k.','markersize',20)
xlim([1 N])
txt = ['$t = 20$'];
T = text(50,8,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients of Temperature with WT $t_{2}$','$z$','$Value$','Coeffcients of Temperature with WT t2',printfigure)

dz = z_2 - z_1;
figure
plot(dz,'k.','markersize',20)
xlim([1 N])
txt = ['$19.5s \to 20s$'];
T = text(50,0.1,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients Change $t_{1}\to t_{2}$','$z$','$Value$','Coeffcients Change with WT t1 to t2',printfigure)

z_sort = sort(abs(dz),'descend');

SW = z_sort(K);
z_k = zeros(size(dz));

for i = 1 : length(dz)
    if abs(dz(i)) >= SW
        z_k(i) = dz(i);
    end
end

figure
plot(z_k,'k.','markersize',20)
xlim([1 N])
txt = ['$19.5s \to 20s$'];
T = text(50,0.1,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('K-term Coeffcients Change with WT','$z$','Value','K-term Coeffcients Change with WT',printfigure)

z_2_k = z_1 + z_k;

y_2_k = THETA * z_2_k;


figure
plot(x,f(:,2001),'k-','linewidth',5)
hold on
plot(x,f(:,2001-50),'c--','linewidth',5)
hold on
plot(xm,y_2_k,'.-','markersize',30)
xlim([0 10])
ylim([0 2])
legend('t = 20','t = 19.5','Signal K-term')
txt = ['$t = 20$'];
T = text(0.8,0.6,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Tempreature Distribution','$x$','$f$','K-term Change Approximation WT',printfigure)




M = 12;
S = zeros(1,M);
S(1:3) = [20,32,45];
k = 4;
while k <= M
    temp = ceil(rand()*N);
    if abs(S - temp) ~= 0
        S(k) = temp;
        k = k + 1;
    end
end
S = sort(S);

Phi = zeros(M,N);
for i = 1 : M
    Phi(i,S(i)) = 1;
end

A = Phi * THETA;

y_sampling = y_2(S);

e = 0.1;
cvx_begin
    variable a(N,1)
    minimize(norm(a,1))
    subject to
        norm(A *(z_1 + a) - y_sampling) <= e
cvx_end


y_2_re = THETA *(z_1 + a);
figure
plot(x,f(:,2001),'k-','linewidth',5)
hold on
plot(x,f(:,2001-50),'c--','linewidth',5)
hold on
plot(xm,y_2_re,'.-','markersize',30)
hold on
plot(xm(S),y_2(S),'r.','Markersize',30)
txt = ['$M = ',num2str(M),'$'];
T = text(0.8,0.8,txt,'FontSize',30);
set(T,'Interpreter','latex')
xlim([0 10])
ylim([0 2])
legend('t = 20','t = 19.5','Signal CS','Measurements')
setplt('Tempreature Distribution','$x$','$f$','CS Change Recovery WT',printfigure)


