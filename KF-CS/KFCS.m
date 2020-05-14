clc
clear
close all

printfigure = 1;

load('Messwerte')

% S = [1,4,7,11];  % benutzte Sensoren
% S = [1,3,5,7,9,11];  % benutzte Sensoren
% S = [1,2,3,4,5,6,7,8,9,10,11];  % benutzte Sensoren
N_time = size(m,1);

number = 6;
S = zeros(N_time,number);
S(:,1) = ceil(rand(N_time,1)*11);
for t = 1 : N_time
    k = 2;
    while k <= number
        temp = ceil(rand()*11);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end

N = 11;
order = 10;
G = 11;      % Grad + 1, anschliesslich 0 Grad.

H = zeros(N, G);  % Eigenfunktionen
lambda = 0 : order;
lambda = lambda * pi / Length;

H(:, 1) = sqrt(1 / Length);
for i = 2 : G
    for n = 1 : N
        H(n, i) = sqrt(2 / Length) * cos(lambda(i) * p(n));
    end
end

A = zeros(G);
for i = 1 : G
    A(i,i) = 1 - step_time * k * lambda(i)^2;
end

Phi = zeros(N_time, number, N);
for t = 1 : N_time
    for i = 1 : number
        Phi(t,i,S(t,i)) = 1;
    end
end

T = zeros(G, N_time);
Ce = zeros(G, G, N_time);
Ce(:,:,1) = eye(G) * 10 ^ 10;
Cv = eye(N) * 1;  % Messunsicherheit
Cw = eye(G) * 1;  % Systemrauschen


for t = 2 : N_time
    h = squeeze(Phi(t,:,:)) * H;
    cv = squeeze(Phi(t,:,:)) * Cv * squeeze(Phi(t,:,:))';
    y = m(t,S(t,:))';
    Tp = A * T(:, t-1);
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * h' / (h * Cp * h' + cv);
    T(:,t) = (eye(size(K,1)) - K * h) * Tp + K * y;
    Ce(:,:,t) = (eye(size(K,1)) - K * h) * Cp;
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
for n = 1 : 400 : N_time
    clf
    plot(x, f(n,:),'k-','LineWidth',5)
    hold on
    plot(x, f_e(n,:),'c-','LineWidth',5)
    hold on
    plot(S(n,:)-1,m(n,S(n,:))','r.','Markersize',40)
    legend('TV Real','TV Estimated','Measure Points')
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
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['KFCS.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

save('KFCS.mat','f_e')


