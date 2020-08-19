clc
clear
close all

printfigure = 1;

load('Messwerte.mat')


y = f(:,end);
x = 0 : dx : 10;
N = length(y);
K= 10;

THETA = DWT(N, 'haar')^-1;

WT = THETA^-1;
z = WT * y;

plot(sort(abs(z),'descend' ),'k.','markersize',20)
xlim([1 N])
txt = ['$t = 20$'];
T = text(800,25,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients of Temperature with WT','$z$','$Value$','Coeffcients of Temperature with WT',printfigure)


z_sort = sort(abs(z),'descend');

SW = z_sort(K);
z_k = zeros(size(z));

for i = 1 : length(z)
    if abs(z(i)) >= SW
        z_k(i) = z(i);
    end
end

figure
plot(sort(abs(z_k),'descend' ),'k.','markersize',20)
xlim([1 N])
txt = ['$t = 20$'];
T = text(800,25,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('K-term Coeffcients of Temperature with WT','$z$','$Value$','K-term Coeffcients of Temperature with WT',printfigure)

y_k = THETA * z_k;

figure
plot(x,y,'k-','linewidth',5)
hold on
plot(x,y_k,'c-','linewidth',5)
xlim([0 10])
ylim([0 2])
legend('Signal','Signal K-term')
txt = ['$t = 20$'];
T = text(0.8,0.6,txt,'FontSize',30);
set(T,'Interpreter','latex')
txt = ['$K = 10$'];
T = text(0.8,0.4,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Tempreature Distribution','$x$','$f$','K-term Approximation WT',printfigure)

M = 24;
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

y_sampling = y(S);

e = 0.03;
cvx_begin
    variable a(N,1)
    minimize(norm(a,1))
    subject to
        norm(A * a - y_sampling) <= e
cvx_end

y_re = THETA * a;

figure
plot(x,y,'k-','linewidth',5)
hold on
plot(x,y_re,'c-','linewidth',5)
hold on
plot(x(S),y_sampling,'r.','markersize',30)
legend('Signal','Signal Estimated','Measure Points')
txt = ['$t = 20$'];
T = text(0.8,0.6,txt,'FontSize',30);
set(T,'Interpreter','latex')
txt = ['$M = ',num2str(M),'$'];
T = text(0.8,0.8,txt,'FontSize',30);
set(T,'Interpreter','latex')
xlim([0 10])
ylim([0 2])
setplt('Recovered Signal by CS','$x$','$f$','Recovered Signal by CS WT',printfigure)
