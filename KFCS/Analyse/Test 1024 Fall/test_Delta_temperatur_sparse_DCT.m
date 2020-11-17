clc
clear
close all

printfigure = 1;

load('Messwerte.mat')

y_1 = f(:,2001-50);
y_2 = f(:,2001);
x = 0 : dx : 10;
N = length(y_1);
K = 10;

THETA = idct(eye(N));

DCT = THETA^-1;
z_1 = DCT * y_1;
z_2 = DCT * y_2;

plot(z_1,'k.','markersize',20)
xlim([1 N])
txt = ['$t = 19.5$'];
T = text(800,25,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients of Temperature with DCT $t_{1}$','$z$','$Value$','Coeffcients of Temperature with DCT t1',printfigure)

figure
plot(z_2,'k.','markersize',20)
xlim([1 N])
txt = ['$t = 20$'];
T = text(800,25,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients of Temperature with DCT $t_{2}$','$z$','$Value$','Coeffcients of Temperature with DCT t2',printfigure)

dz = z_2 - z_1;
figure
plot(dz,'k.','markersize',20)
xlim([1 N])
txt = ['$19.5s \to 20s$'];
T = text(800,0.2,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients Change $t_{1}\to t_{2}$','$z$','$Value$','Coeffcients Change with DCT t1 to t2',printfigure)

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
txt = ['$t = 20$'];
T = text(800,25,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('K-term Coeffcients Change with DCT','$z$','$Value$','K-term Coeffcients Change with DCT',printfigure)

z_2_k = z_1 + z_k;

y_2_k = THETA * z_2_k;


figure
plot(x,y_2,'k-','linewidth',5)
hold on
plot(x,y_1,'b-','linewidth',5)
hold on
plot(x,y_2_k,'c-','linewidth',5)
xlim([0 10])
ylim([0 2])
legend('t = 20','t = 19.5','Signal K-Terme')
txt = ['$t = 20$'];
T = text(0.8,0.8,txt,'FontSize',60);
set(T,'Interpreter','latex')
txt = ['$K = 10$'];
T = text(0.8,0.4,txt,'FontSize',60);
set(T,'Interpreter','latex')
setplt('Temperaturverteilung','$x$','$f$','K-term Change Approximation DCT',printfigure)




M = 12;
S = zeros(1,M);
S(1) = ceil(rand()*N);
k = 2;
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

e = 0.01;
cvx_begin
    variable a(N,1)
    minimize(norm(a,1))
    subject to
        norm(A *(z_1 + a) - y_sampling) <= e
cvx_end


y_2_re = THETA *(z_1 + a);
figure
plot(x,y_2,'k-','linewidth',5)
hold on
plot(x,y_1,'b-','linewidth',5)
hold on
plot(x,y_2_re,'c--','linewidth',5)
hold on
plot(x(S),y_2(S),'r.','Markersize',40)
xlim([0 10])
ylim([0 2])
txt = ['$M = 24$'];
T = text(0.8,0.8,txt,'FontSize',60);
set(T,'Interpreter','latex')
legend('t = 20','t = 19.5','Signal aus CS','Messungen')
setplt('Rekonstruiertes Signal aus CS','$x$','$f$','CS Change Recovery DCT',printfigure)

close all
