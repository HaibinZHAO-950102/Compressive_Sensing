clc
clear
close all

printfigure = 1;

load('Messwerte.mat')


y = f(p_index,end);
x = 0 : dx : 10;
xm = x(p_index);
N = length(y);
K= 10;

DCT = idct(eye(N))^-1;
THETA = DCT^-1;

WT = THETA^-1;
z = WT * y;

plot(sort(abs(z),'descend' ),'k.','markersize',20)
xlim([1 N])
txt = ['$t = 20$'];
T = text(800,25,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Coeffcients of Temperature with DCT','$z$','Value','Coeffcients of Temperature with DCT',printfigure)


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
setplt('K-term Coeffcients of Temperature with DCT','$z$','Value','K-term Coeffcients of Temperature with DCT',printfigure)

y_k = THETA * z_k;

figure
plot(x,f(:,end),'k-','linewidth',3)
hold on
plot(xm,y_k,'.-','markersize',30)
xlim([0 10])
ylim([0 2])
legend('Signal','Signal K-term')
txt = ['$t = 20$'];
T = text(0.8,0.6,txt,'FontSize',30);
set(T,'Interpreter','latex')
setplt('Tempreature Distribution','$x$','$f$','K-term Approximation DCT',printfigure)

M = 20;
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
plot(x,f(:,end),'k-','linewidth',3)
hold on
plot(xm,y_re,'.-','markersize',30)
hold on
plot(xm(S),y_sampling,'r.','markersize',20)
legend('Signal','Signal Estimated','Measure Points')
txt = ['$t = 20$'];
T = text(0.8,0.6,txt,'FontSize',30);
set(T,'Interpreter','latex')
txt = ['$M = ',num2str(M),'$'];
T = text(0.8,0.8,txt,'FontSize',30);
set(T,'Interpreter','latex')
xlim([0 10])
ylim([0 2])
setplt('Recovered Signal by CS','$x$','$f$','Recovered Signal by CS DCT',printfigure)
