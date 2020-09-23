
clc
clear
close all

printfigure = 1;

load('Messwerte')

Dt = 0.1;
m = m(:,1:Dt/dt:end);

N_time = size(m,2);


number = 12;
S = zeros(N_time,number);
S(:,1) = ceil(rand(N_time,1)*64);
for t = 1 : N_time
    kk = 2;
    while kk <= number
        temp = ceil(rand()*64);
        if abs(S(t,:) - temp) ~= 0
            S(t,kk) = temp;
            kk = kk + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end

N = 64;
order = 50;
G = order+1;      % Grad + 1, anschliesslich 0 Grad.

H = zeros(N, G);  % Eigenfunktionen
lambda = 0 : order;
lambda = lambda * pi / Length;

H(:, 1) = sqrt(1 / Length);
for i = 2 : G
    for n = 1 : N
        H(n, i) = sqrt(2 / Length) * cos(lambda(i) * p(n));
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


Phi = zeros(N_time, number, N);
for t = 1 : N_time
    for i = 1 : number
        Phi(t,i,S(t,i)) = 1;
    end
end

T = zeros(G, N_time);
T(:,1) = (H'*H)^-1 * H' * m(:,1);
Ce = zeros(G, G, N_time);
Ce(:,:,1) = eye(G) * 10 ^ 10;
Cv = eye(N) * 0.01;  % Messunsicherheit
Cw = eye(G) * 1;  % Systemrauschen


for t = 2 : N_time
    h = squeeze(Phi(t,:,:)) * H;
    cv = squeeze(Phi(t,:,:)) * Cv * squeeze(Phi(t,:,:))';
    y = m(S(t,:),t);
    Tp = A * T(:, t-1);
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * h' / (h * Cp * h' + cv);
    T(:,t) = (eye(size(K,1)) - K * h) * Tp + K * y;
    Ce(:,:,t) = (eye(size(K,1)) - K * h) * Cp;
end

x = 0 : dx : Length;
N_length = length(x);
f_e = zeros(N_length,N_time);  % Temperaturmatrix
phi = zeros(N_length,G);  % Eigenfunktionen
phi(:,1) = sqrt(1 / Length);
for i = 2 : G
    for n = 1 : N_length
        phi(n,i) = sqrt(2 / Length) * cos(lambda(i) * x(n));
    end
end
for n = 1 : N_time
    f_e(:,n) = phi * T(:,n);
end

figure
for n = 1 : 0.5/Dt : N_time
    clf
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    plot(x, f(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e(:,n),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(n,:))),m(S(n,:),n)','r.','Markersize',40)
    legend('Signal','Signal Estimated','Measure Points')
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
    name = ['random_simpling_KF_modal.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% f_random_simpling_KF_modal = f_e;
% save('f_random_simpling_KF_modal.mat','f_random_simpling_KF_modal')


