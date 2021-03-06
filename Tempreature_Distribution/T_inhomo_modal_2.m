clc
clear
close all

printfigure = 1;

Length = 10;  % Stablaenge
Time = 20;   % Zetiraum

dx = Length / 1023;
dt = 0.01;
x = 0 : dx : Length;
t = 0 : dt : Time;
nx = length(x);
nt = length(t);

k = 0.1; % Waermeleitfaehigkeit in cm^2/s

N = 50;  % Grad
lambda = 0 : pi / Length : N * pi / Length;
f = zeros(nt, nx);  % Temperaturmatrix
f(1,:) = sin(x / Length * 2 * pi) + 1;
% f(1,:) = floor(x * 5) / 5;
phi = zeros(N + 1, nx);  % Eigenfunktionen
T = zeros(N + 1, nt);  % Gewichtung
u = zeros(nt, nx);  % Anregung
U = zeros(N + 1, nt);  % Anregungsgewichtung

plot(x,f(1,:),'LineWidth',5)
setplt('Anfangsbedingung','$x$','$f$','TV_inhomo_modal_inital_condition',0)

phi(1,:) = sqrt(1 / Length);
for i = 2 : N + 1
    for n = 1 : nx
        phi(i, n) = sqrt(2 / Length) * cos(lambda(i) * x(n));
    end
end

figure
for i = 1 : N + 1
    plot(0:dx:Length,phi(i,:))
    hold on
end
txt = ['$G = ',num2str(N),'$'];
TEXT = text(8,0.4,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunktionen','$x$','','TV_inhomo_modal_Eigenfunctions',0)



u(:,round(3/dx)+1) = 0.1 * sin(t - pi / 4);
u(:,round(5/dx)+1) = -0.2 * sin(t);
u(:,round(7/dx)+1) = 0.01 * t;

for i = 1 : N + 1
    U(i,:) = u(:,round(3/dx)+1) * phi(i,round(3/dx)+1) + u(:,round(5/dx)+1) * phi(i,round(5/dx)+1) + u(:,round(7/dx)+1) * phi(i,round(7/dx)+1);
end

for i = 1 : N + 1
    T(i, 1) = 0;
    for n = 1 : nx
        T(i, 1) = T(i, 1) + f(1, n) * phi(i, n) * dx;
    end
end

for i = 1 : N + 1
    for n = 2 : nt
        T(i, n) = (1 - dt * k * lambda(i)^2) * T(i, n - 1) + dt * U(i, n);
    end
end

f = zeros(nt, nx);  % Temperaturmatrix
for n = 1 : nt
    for i = 1 : N + 1
        f(n,:) = f(n,:) + T(i,n) * phi(i,:);
    end
end

figure
for n = 1 : 0.5/dt : nt
    plot(x, f(n,:),'LineWidth',5)
    ylim([-0.5 2.5])
    setplt('','$x$','$f$','Temperature Distribution',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*dt),'$'];
    TEXT = text(8,1.4,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$G = ',num2str(N),'$'];
    TEXT = text(8,1.8,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'TV_inhomo_modal_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'TV_inhomo_modal_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

figure
[X, Y] = meshgrid(x, t);
mesh(X,Y,f)
setmesh('','$x$','$t$','$f$','TV_inhomo_modal_Ttx_2',printfigure)


% Measurement
% positions
p_index = round(linspace(1,nx,64)); % 64 sensors from nx points
p = x(p_index);                     % positions of 51 sensors from 0 - Length
m = f(:,p_index)';                  % their measurements

% f_modal_inhomo_2 = f;
% x_modal_inhomo_2 = x;
% save('f_modal_inhomo_2.mat','x_modal_inhomo_2','f_modal_inhomo_2')
% 
% f = f';
% save('Messwerte.mat','f','k','Length','dt','dx','p','m','p_index')

close all
