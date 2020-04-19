clc
clear
close all

printfigure = 0;

load('Messwerte')

% S = [1,4,7,11];  % benutzte Sensoren
% S = [1,3,5,7,9,11];  % benutzte Sensoren
S = [1,2,3,4,5,6,7,8,9,10,11];  % benutzte Sensoren

N = length(S);  % Anzahl der Messungen
order = N - 1;  % Grad
G = order + 1;      % Grad + 1, anschliesslich 0 Grad.

N_time = size(m,1);
H = zeros(N, G);  % Eigenfunktionen
lambda = 0 : order;
lambda = lambda * pi / Length;

H(:, 1) = sqrt(1 / Length);
for i = 2 : G
    for n = 1 : N
        H(n, i) = sqrt(2 / Length) * cos(lambda(i) * p(S(n)));
    end
end

A = zeros(G);
for i = 1 : G
    A(i,i) = 1 - step_time * k * lambda(i)^2;
end

T = zeros(G, N_time);
Ce = zeros(G, G, N_time);
Ce(:,:,1) = eye(G) * 1000000000000;
Cv = eye(N) * 1;  % Messunsicherheit
Cw = eye(G) * 1;  % Systemrauschen
for t = 2 : N_time
    Tp = A * T(:, t-1);
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,t) = (eye(size(K,1)) - K * H) * Tp + K * m(t,S)';
    Ce(:,:,t) = (eye(size(K,1)) - K * H) * Cp;
end

step_length = 0.01;
x = 0 : step_length : Length;
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
for n = 1 : 40 : N_time
    clf
    plot(x, f(n,:),'k-','LineWidth',5)
    hold on
    plot(x, f_e(n,:),'c-','LineWidth',5)
    hold on
    plot(p(S),f(n,p(S)/step_length+1),'r.','Markersize',40)
    legend('TV Real','TV Estimated','Measure Points')
    xlim([0 10])
    ylim([0 2])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = ',num2str(N),'$'];
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




