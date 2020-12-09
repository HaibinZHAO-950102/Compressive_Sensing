clc
clear
close all

printfigure = 1;

Length_x = 10;  
Length_y = 20;  
Time = 200;   % Zetiraum
step_length_x = 0.1;
step_length_y = 0.2;
step_time = 1;
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
        f(1,i,j) = sin(x(i) / Length_x * 2 * pi) * cos(y(j) / Length_y * 2 * pi);
    end
end

[Y,X] = meshgrid(y,x);
mesh(X,Y,squeeze(f(1,:,:)));
caxis([-1 1])
pbaspect([1 Length_y/Length_x 0.5])
setmesh('Anfangsbedingung','$x$','$y$','$f$','T_2D_homo_modal_inital_condition',printfigure)

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
setplt('Eigenfunktionen $\varphi$','$x$','$value$','T_2D_homo_modal_Eigenfunctions_phi',printfigure)

figure
for i = 1 : N + 1
    plot(y,psi(i,:))
    hold on
end
txt = ['$G = ',num2str(N),'$'];
TEXT = text(16,0.3,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunktionen $\psi$','$y$','$value$','T_2D_homo_modal_Eigenfunctions_psi',printfigure)

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
            T(i, j, n) = (1 - step_time * k * (lambda(i)^2 + sigma(j)^2))^(n - 1) * T(i, j, 1);
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
for n = 1 : 5 : N_time
    mesh(X,Y,squeeze(f(n,:,:)))
    zlim([-1.5 1.5])
    xticks(0:5:10)
    yticks(0:5:20)
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    pbaspect([1 Length_y/Length_x 0.5])
    setmesh('','$x$','$y$','$f$','T_2D_homo_modal',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    TEXT = text(8,0,0.5,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$G = ',num2str(N),'$'];
    TEXT = text(8,0,1.8,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    caxis([-1 1])
    
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'T_2D_homo_modal.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'T_2D_homo_modal.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

figure
for n = 1 : 6
    a = round((N_time - 1) / 5 * (n - 1)) + 1;
    mesh(X,Y,squeeze(f(a,:,:)))
    zlim([-1.5 1.5])
    pbaspect([1 Length_y/Length_x 0.5])
    caxis([-1 1])
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(a)),'$'];
    TEXT = text(8,0,0.5,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$G = ',num2str(N),'$'];
    TEXT = text(8,0,1.8,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    drawnow
    name = ['T_2D_homo_modal_shot_',num2str(n)];
    setmesh('','$x$','$y$','$f$',name,printfigure)
end

close all