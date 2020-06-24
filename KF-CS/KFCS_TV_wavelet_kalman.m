% wavelet + kalman

clc
clear
close all

printfigure = 0;

load('kfcs_tv_wavelet_6_l2');
load('Messwerte')

x = 0 : 0.01 : 10-0.01;
nx = length(x);
t = 0 : 0.01 : 20;
nt = length(t);

index = 1:10:nx;

m = f_e_wt(index,:);
xm = x(index);

figure
for n = 1 : 40 : nt
    clf
    plot(x, f_e_wt(:,n),'k-','LineWidth',3)
    hold on
    plot(xm,m(:,n)','r.','Markersize',40)
    legend('WT Estimated Signal','Measure Poins')
    xlim([0 10])
    ylim([-0.5 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',30);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*0.01),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = 6$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_wavelet_selected_points.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



lambda = 0 : 50;
lambda = lambda * pi / 10;

Psi = zeros(nx, 51);  % Eigenfunktionen

Psi(:, 1) = sqrt(1 / 10);
for i = 2 : 51
    for n = 1 : nx
        Psi(n, i) = sqrt(2 / 10) * cos(lambda(i) * x(n));
    end
end

Phi = zeros(length(index),nx);
for i = 1 : length(index)
    Phi(i,index(i)) = 1;
end

H = Phi * Psi;

k = 0.1;

A = zeros(51);
for i = 1 : 51
    A(i,i) = 1 - 0.01 * k * lambda(i)^2;
end

T = zeros(51, nt);
Ce = zeros(51, 51, nt);
Ce(:,:,1) = eye(51) * 10 ^ 10;
Cv = eye(length(index)) * 10;  % Messunsicherheit
Cw = eye(51) * 1;  % Systemrauschen

for n = 2 : nt
    Tp = A * T(:, n-1);
    Cp = A * Ce(:,:, n-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,n) = (eye(size(K,1)) - K * H) * Tp + K * m(:,n);
    Ce(:,:,n) = (eye(size(K,1)) - K * H) * Cp;
end

f_e_wt_kf = zeros(size(f_e_wt));
f_e_wt_kf(:,1) = f_e_wt(:,1);
for i = 2 : nt
    f_e_wt_kf(:,i) = Psi * T(:,i);
end

figure
for n = 1 : 40 : nt
    clf
    plot(x, f(n,1:end-1),'k-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf(:,n),'c-','LineWidth',3)
    hold on
    plot(x((S(n,:)-1)*10+1),f(n,(S(n,:)-1)*10+1)','r.','Markersize',40)
    legend('Signal','Signal Estimated','Measure Poins')
    xlim([0 10])
    ylim([-0.5 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*0.01),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = 6$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_wavelet_KF.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

save('f_e_wt_kf.mat','f_e_wt_kf')







