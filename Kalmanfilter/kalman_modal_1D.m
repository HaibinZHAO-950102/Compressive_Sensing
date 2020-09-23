clc
clear
close all

printfigure = 0;

load('Messwerte')

N = 36;  % Anzahl der Messungen
Dt = 0.1; % time_step

S = round(linspace(1,64,N));  % benutzte Sensoren
m = m(S,1:Dt/dt:end);

order = N - 1;  % Grad
G = order + 1;      % Grad + 1, anschliesslich 0 Grad.



N_time = size(m,2);
H = zeros(N, G);  % Eigenfunktionen
lambda = 0 : order;
lambda = lambda * pi / Length;

H(:, 1) = sqrt(1 / Length);
for i = 2 : G
    for n = 1 : N
        H(n, i) = sqrt(2 / Length) * cos(lambda(i) * p(S(n)));
    end
end

Dt_max = 2 / k / lambda(end)^2;
if Dt > Dt_max
    ['Dt should be under ',num2str(Dt_max)]
end

A = zeros(G);
for i = 1 : G
    A(i,i) = 1 - Dt * k * lambda(i)^2;
end



T = zeros(G, N_time);
T(:,1) = (H'*H)^-1 * H' * m(:,1);


Ce = zeros(G, G, N_time);
Ce(:,:,1) = eye(G) * 10 ^ 10;
Cv = eye(N) * 1;  % Messunsicherheit
Cw = eye(G) * 1;  % Systemrauschen
for t = 2 : N_time
    Tp = A * T(:, t-1);
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,t) = (eye(size(K,1)) - K * H) * Tp + K * m(:,t);
    Ce(:,:,t) = (eye(size(K,1)) - K * H) * Cp;
end

x = 0 : dx : Length;
N_length = length(x);
f_e = zeros(N_time, N_length);  % Temperaturmatrix
phi = zeros(G, N_length);  % Eigenfunktionen
phi(1,:) = sqrt(1 / Length);
for i = 2 : G
    for n = 1 : N_length
        phi(i, n) = sqrt(2 / Length) * cos(lambda(i) * x(n));
    end
end
for n = 1 : N_time
    for i = 1 : G
        f_e(n,:) = f_e(n,:) + T(i,n) * phi(i,:);
    end
end

figure
for n = 1 : 0.5/Dt : N_time
    clf
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    plot(x, f(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e(n,:),'c-','LineWidth',5)
    hold on
    plot(p(S),m(:,n),'r.','Markersize',40)
    legend('Signal','Signal Estimated','Measurements')
    xlim([0 10])
    ylim([-0.5 2.5])
    setplt('Temperature Distribution','$x$','$f$','Temperature Distribution',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(N),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['Kalman_modal_',num2str(N),'.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% f_e_kalman_modal_1D_36 = f_e';
% save('f_e_kalman_modal_1D_36.mat','f_e_kalman_modal_1D_36')

