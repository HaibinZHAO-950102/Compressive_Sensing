% dct + kalman integrated 2

clc
clear
close all

printfigure = 0;

load('Messwerte')

x = 0 : dx : 10;
nx = length(x);

Dt = 0.1;
m = m(:,1:Dt/dt:end);

nt = size(m,2);

% random select of sensors
number = 12;
S = zeros(nt,number);
S(:,1) = ceil(rand(nt,1)*64);
for t = 1 : nt
    k = 2;
    while k <= number
        temp = ceil(rand()*64);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end



% CS Basis

DCT = idct(eye(64))^-1;
THETA = DCT^-1;


Phi_cs = zeros(nt, number, 64);
for t = 1 : nt
    for i = 1 : number
        Phi_cs(t,i,S(t,i)) = 1;
    end
end

z = zeros(64, nt);
z(:,1) = DCT * m(:,1);


% Kalman-modell
G = 20; % order

lambda = 0 : G;
lambda = lambda * pi / 10;

Psi = zeros(nx, G+1);  % Eigenfunktionen

Psi(:, 1) = sqrt(1 / 10);
for i = 2 : G+1
    for n = 1 : nx
        Psi(n, i) = sqrt(2 / 10) * cos(lambda(i) * x(n));
    end
end

Phi_KF = zeros(length(p_index),nx);
for i = 1 : length(p_index)
    Phi_KF(i,p_index(i)) = 1;
end

H = Phi_KF * Psi;

k = 0.1;

A = zeros(G+1);
for i = 1 : G+1
    A(i,i) = 1 - Dt * k * lambda(i)^2;
end

T = zeros(G+1, nt);
Ce = zeros(G+1, G+1, nt);
Ce(:,:,1) = eye(G+1) * 10 ^ 10;
Cv = eye(length(p_index)) * 20;  % Messunsicherheit
Cw = eye(G+1) * 5;  % Systemrauschen

f_e_dct_kf_integrated_2 = zeros(nx,nt);
f_e_dct_kf_integrated_2(:,1) = f(:,1);

dz = zeros(64,nt);

T(:,1) = (Psi'*Psi)^-1*Psi'*f_e_dct_kf_integrated_2(:,1);

% KFCS
for t = 2 : nt
    t
    % CS estimation
    h = squeeze(Phi_cs(t,:,:)) * THETA;
    y = m(S(t,:),t);
    
    Tp = A * T(:, t-1);
    zp =  DCT * Phi_KF * Psi * Tp;
    
    e = 0;
    cvx_begin quiet;
        variable a(64,1);
        minimize(norm(a,1));
        subject to;
            norm(h * (zp + a) - y) <= e;
    cvx_end;
    z(:, t) = zp + a;
    dz(:,t-1) = a;
    yp =  THETA * z(:, t);
   
    Cv = eye(length(p_index)) * 20; 
    for i = 1 : size(S,2)
        Cv(S(t,i),S(t,i)) = 0.1;
    end

    
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,t) = (eye(size(K,1)) - K * H) * Tp + K * yp;
    Ce(:,:,t) = (eye(size(K,1)) - K * H) * Cp;
end

for i = 1 : nt
    f_e_dct_kf_integrated_2(:,i) = Psi * T(:,i);
end

f = f(:,1:Dt/dt:end);
figure
for n = 1 : 0.5/Dt : nt
    clf
    plot(x, f(:,n),'k-','LineWidth',3)
    hold on
    plot(x, f_e_dct_kf_integrated_2(:,n),'c-','LineWidth',3)
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
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_wavelet_KF_modal_integrated_2.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

f_e_dct_kf_modal_integrated_2 = f_e_dct_kf_integrated_2;
save('f_e_dct_kf_modal_integrated_2.mat','f_e_dct_kf_modal_integrated_2')


