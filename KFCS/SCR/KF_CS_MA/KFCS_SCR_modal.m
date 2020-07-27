% mt + kalman

clc
clear
close all

printfigure = 1;

load('CS_SCR_L1');
load('Messwerte')

x = 0 : dx : 10;
nx = length(x);
Dt = 0.1;
t = 0 : Dt : 20;
nt = length(t);


m = yp;
xm = x(p_index);

G = 50;
lambda = 0 : G;
lambda = lambda * pi / 10;

Psi = zeros(nx, G+1);  % Eigenfunktionen

Psi(:, 1) = sqrt(1 / 10);
for i = 2 : G+1
    for n = 1 : nx
        Psi(n, i) = sqrt(2 / 10) * cos(lambda(i) * x(n));
    end
end

Phi = zeros(length(p_index),nx);
for i = 1 : length(p_index)
    Phi(i,p_index(i)) = 1;
end

H = Phi * Psi;

k = 0.1;

A = zeros(G+1);
for i = 1 : G+1
    A(i,i) = 1 - Dt * k * lambda(i)^2;
end

T = zeros(G+1, nt);
Ce = zeros(G+1, G+1, nt);
Ce(:,:,1) = eye(G+1) * 10 ^ 10;
Cv = eye(length(p_index)) * 0.05;  % Messunsicherheit
Cw = eye(G+1) * 2;  % Systemrauschen

for n = 2 : nt
    % dynamic gewichtung
    Cv = eye(length(p_index)) * 2; 
    for i = 1 : size(S,2)
        Cv(S(n,i),S(n,i)) = 0.05;
    end
    
    Tp = A * T(:, n-1);
    Cp = A * Ce(:,:, n-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,n) = (eye(size(K,1)) - K * H) * Tp + K * m(:,n);
    Ce(:,:,n) = (eye(size(K,1)) - K * H) * Cp;
end

f_e_scr_kf = zeros(nx,nt);
for i = 2 : nt
    f_e_scr_kf(:,i) = Psi * T(:,i);
end

figure
for n = 1 : 0.5/Dt : nt
    clf
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    plot(x, f(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e_scr_kf(:,n),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(n,:))),m(S(n,:),n),'r.','Markersize',40)
    legend('Signal','Signal Estimated','Measure Poins')
    xlim([0 10])
    ylim([-0.5 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$f$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(size(S,2)),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_SCR_KF_modal.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

f_e_SCR_kf_modal = f_e_scr_kf;
save('f_e_SCR_kf_modal.mat','f_e_SCR_kf_modal')







