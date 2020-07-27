clc
clear
close all

printfigure = 1;

Length_x = 10;  
Length_y = 20;  
Time = 20;
step_length_x = 0.1;
step_length_y = 0.2;
step_time = 0.4;
x = 0 : step_length_x : Length_x;
y = 0 : step_length_y : Length_y;
t = 0 : step_time : Time;
N_length_x = length(x);
N_length_y = length(y);
N_time = length(t);

k = 0.1;  % Waermeleitfaehigkeit in cm^2/s

N = 10;  % Grad
lambda = 0 : pi / Length_x : N * pi / Length_x;
sigma = 0 : pi / Length_y : N * pi / Length_y;

f = zeros(N_time, N_length_x, N_length_y);  % Temperaturmatrix
for i = 1 : N_length_x
    for j = 1 : N_length_y
        f(1,i,j) = sin(x(i) / Length_x * 2 * pi) * cos(y(j) / Length_y * 2 * pi) + 1;
    end
end

[Y,X] = meshgrid(y,x);
mesh(X,Y,squeeze(f(1,:,:)));
caxis([-0.5 2.5])
pbaspect([1 Length_y/Length_x 0.5])
setmesh('Initial Condition','$x$','$y$','$f$','T_2D_inhomo_modal_inital_condition',0)

phi = zeros(N + 1, N_length_x);  % Eigenfunktionen
psi = zeros(N + 1, N_length_y);  % Eigenfunktionen
T = zeros(N + 1, N + 1, N_time);  % Gewichtung

phi(1,:) = sqrt(1 / Length_x);
for i = 2 : N + 1
    for n = 1 : N_length_x
        phi(i, n) = sqrt(2 / Length_x) * cos(lambda(i) * x(n));
    end
end
psi(1,:) = sqrt(1 / Length_y);
for i = 2 : N + 1
    for n = 1 : N_length_y
        psi(i, n) = sqrt(2 / Length_y) * cos(sigma(i) * y(n));
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

u = zeros(N_time, N_length_x, N_length_y);  % Anregung
U = zeros(N + 1, N + 1, N_time);

u(:,31,16) = 1 * sin(t - pi / 4);
u(:,51,51) = -2 * sin(t);
u(:,71,66) = 0.1 * t;
u(:,81,91) = 0.3 * t;

for m = 1 : N + 1
    for n = 1 : N + 1
        for i = 1 : N_time
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
        for n = 1 : N_length_x
            for m = 1 : N_length_y
                T(i, j, 1) = T(i, j, 1) + f(1,n,m) * phi(i, n) * psi(j, m) * step_length_x * step_length_y;
            end
        end
    end
end

for i = 1 : N + 1
    for j = 1 : N + 1
        for n = 2 : N_time
            T(i, j, n) = (1 - step_time * k * (lambda(i)^2 + sigma(j)^2)) * T(i, j, 1) + step_time * U(i, j, n);
        end
    end
end

f = zeros(N_time, N_length_x, N_length_y);  % Temperaturmatrix
for n = 1 : N_time
    for i = 1 : N + 1
        for j = 1 : N + 1
            for x_direction = 1 : N_length_x
                for y_direction = 1 : N_length_y
                    f(n,x_direction,y_direction) = f(n,x_direction,y_direction) + T(i,j,n) * phi(i,x_direction) * psi(j,y_direction);
                end
            end
        end
    end
end

figure
for n = 1 : 1 : N_time
    mesh(X,Y,squeeze(f(n,:,:)))
    zlim([-0.5 2.5])
    caxis([-0.5 2.5])
    pbaspect([1 Length_y/Length_x 0.5])
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
    a = round((N_time - 1) / 5 * (n - 1)) + 1;
    mesh(X,Y,squeeze(f(a,:,:)))
    zlim([-0.5 2.5])
    caxis([-0.5 2.5])
    pbaspect([1 Length_y/Length_x 0.5])
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

