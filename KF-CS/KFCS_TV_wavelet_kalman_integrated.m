% wavelet + kalman integrated

clc
clear
close all

printfigure = 0;

load('Messwerte')

x = 0 : 0.01 : 10 - 0.01;
nx = length(x);
p = 0 : 0.1 : 10 - 0.1;
m = f(:,1:10:1000);

nt = size(m,1);

% random select of sensors
number = 6;
S = zeros(nt,number);
S(:,1) = ceil(rand(nt,1)*100);
for t = 1 : nt
    k = 2;
    while k <= number
        temp = ceil(rand()*100);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end


% CS Basis

WT = DWT(1000,'haar');
THETA = WT^-1;


Phi_cs = zeros(nt, number, nx);
for t = 1 : nt
    for i = 1 : number
        Phi_cs(t,i,(S(t,i)-1)*10+1) = 1;
    end
end

z = zeros(nx, nt);
z(:,1) = WT * f(1,1:end-1)';


% Kalman-modell

lambda = 0 : 50;
lambda = lambda * pi / 10;

Psi = zeros(nx, 51);  % Eigenfunktionen

Psi(:, 1) = sqrt(1 / 10);
for i = 2 : 51
    for n = 1 : nx
        Psi(n, i) = sqrt(2 / 10) * cos(lambda(i) * x(n));
    end
end

index = 1 : 10 : nx;
Phi_KF = zeros(length(index),nx);
for i = 1 : length(index)
    Phi_KF(i,index(i)) = 1;
end

H = Phi_KF * Psi;

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

f_e_wt_kf_integrated = zeros(nx,nt);
f_e_wt_kf_integrated(:,1) = THETA * z(:,1);

dz = zeros(1000,nt);

% KFCS
for t = 2 : nt
    % CS estimation
    h = squeeze(Phi_cs(t,:,:)) * THETA;
    y = m(t,S(t,:))';
    
    e = 0;
    cvx_begin;
        variable a(nx,1);
        minimize(norm(a,2));
        subject to;
            norm(h * (z(:,t-1) + a) - y) <= e;
    cvx_end;
    z(:, t) = z(:,t-1) + a;
    dz(:,t-1) = a;
    
    Tp = A * T(:, t-1);
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,t) = (eye(size(K,1)) - K * H) * Tp + K * m(t,:)';
    Ce(:,:,t) = (eye(size(K,1)) - K * H) * Cp;

    f_e_wt_kf_integrated(:,t) = Psi * T(:,t);
    
    z(:, t) = WT * f_e_wt_kf_integrated(:,t);
end

figure
for n = 1 : 40 : nt
    clf
    plot(x, f(n,1:end-1),'k-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf_integrated(:,n),'c-','LineWidth',3)
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
    name = ['CS_wavelet_KF_integrated.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% save('f_e_wt_kf_integrated.mat','f_e_wt_kf_integrated')



