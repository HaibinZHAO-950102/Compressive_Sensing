clc
clear
close all

printfigure = 0;

Lx = 10;  
Ly = 20;  
Time = 20;
dx = Lx / 127;
dy = Ly / 127;
dt = 0.25;
x = 0 : dx : Lx;
y = 0 : dy : Ly;
t = 0 : dt : Time;
nx = length(x);
ny = length(y);
nt = length(t);

k = 0.1;  % Waermeleitfaehigkeit in cm^2/s

N = 10;  % Grad
lambda = 0 : pi / Lx : N * pi / Lx;
sigma = 0 : pi / Ly : N * pi / Ly;

f = zeros(nt, nx, ny);  % Temperaturmatrix
for i = 1 : nx
    for j = 1 : ny
        f(1,i,j) = sin(x(i) / Lx * 2 * pi) * cos(y(j) / Ly * 2 * pi) + 1;
    end
end

[Y,X] = meshgrid(y,x);
mesh(X,Y,squeeze(f(1,:,:)));
caxis([-0.5 2.5])
pbaspect([1 Ly/Lx 0.5])
setmesh('Initial Condition','$x$','$y$','$f$','T_2D_inhomo_modal_inital_condition',0)

phi = zeros(N + 1, nx);  % Eigenfunktionen
psi = zeros(N + 1, ny);  % Eigenfunktionen
T = zeros(N + 1, N + 1, nt);  % Gewichtung

phi(1,:) = sqrt(1 / Lx);
for i = 2 : N + 1
    for n = 1 : nx
        phi(i, n) = sqrt(2 / Lx) * cos(lambda(i) * x(n));
    end
end
psi(1,:) = sqrt(1 / Ly);
for i = 2 : N + 1
    for n = 1 : ny
        psi(i, n) = sqrt(2 / Ly) * cos(sigma(i) * y(n));
    end
end

figure
for i = 1 : N + 1
    plot(x,phi(i,:))
    hold on
end
txt = ['$G = ',num2str(N),'$'];
TEXT = text(8,0.4,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunctions $\varphi$','$x$','$value$','T_2D_inhomo_modal_Eigenfunctions_phi',0)

figure
for i = 1 : N + 1
    plot(y,psi(i,:))
    hold on
end
txt = ['$G = ',num2str(N),'$'];
TEXT = text(16,0.3,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunctions $\psi$','$y$','$value$','T_2D_inhomo_modal_Eigenfunctions_psi',0)

u = zeros(nt, nx, ny);  % Anregung
U = zeros(N + 1, N + 1, nt);

u(:,31,16) = 1 * sin(t - pi / 4);
u(:,51,51) = -2 * sin(t);
u(:,71,66) = 0.1 * t;
u(:,81,91) = 0.3 * t;

for m = 1 : N + 1
    for n = 1 : N + 1
        for i = 1 : nt
            U(m,n,i) = u(i,31,16) * phi(m,31) * psi(n,16) +...
                       u(i,51,51) * phi(m,51) * psi(n,51) +...
                       u(i,71,66) * phi(m,66) * psi(n,66) + ...
                       u(i,81,91) * phi(m,81) * psi(n,91);
        end
    end
end

for i = 1 : N + 1
    for j = 1 : N + 1
        T(i, j, 1) = 0;
        for n = 1 : nx
            for m = 1 : ny
                T(i, j, 1) = T(i, j, 1) + f(1,n,m) * phi(i, n) * psi(j, m) * dx * dy;
            end
        end
    end
end

for i = 1 : N + 1
    for j = 1 : N + 1
        for n = 2 : nt
            T(i, j, n) = (1 - dt * k * (lambda(i)^2 + sigma(j)^2)) * T(i, j, 1) + dt * U(i, j, n);
        end
    end
end

f = zeros(nt, nx, ny);  % Temperaturmatrix
for n = 1 : nt
    for i = 1 : N + 1
        for j = 1 : N + 1
            for x_direction = 1 : nx
                for y_direction = 1 : ny
                    f(n,x_direction,y_direction) = f(n,x_direction,y_direction) + T(i,j,n) * phi(i,x_direction) * psi(j,y_direction);
                end
            end
        end
    end
end

figure
for n = 1 : 2 : nt
    mesh(X,Y,squeeze(f(n,:,:)))
    zlim([-0.5 2.5])
    caxis([-0.5 2.5])
    pbaspect([1 Ly/Lx 0.5])
    setmesh('Temperature Distribution','$x$','$y$','$f$','T_2D_inhomo_modal_2',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(n)),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    txt = ['$G = ',num2str(N),'$'];
    TEXT = text(8,0,1.5,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'T_2D_inhomo_modal_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'T_2D_inhomo_modal_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

figure
for n = 1 : 6
    a = round((nt - 1) / 5 * (n - 1)) + 1;
    mesh(X,Y,squeeze(f(a,:,:)))
    zlim([-0.5 2.5])
    caxis([-0.5 2.5])
    pbaspect([1 Ly/Lx 0.5])
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(a)),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    txt = ['$G = ',num2str(N),'$'];
    TEXT = text(8,0,1.5,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    drawnow
    name = ['T_2D_inhomo_modal_2_shot_',num2str(n)];
    setmesh('Temperature Distribution','$x$','$y$','$T$',name,printfigure)
end

