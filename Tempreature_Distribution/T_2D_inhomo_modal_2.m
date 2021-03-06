clc
clear
close all

printfigure = 1;

Lx = 10;  
Ly = 20;  
Time = 20;
dx = Lx / 127;
dy = Ly / 127;
dt = 0.1;
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

f = zeros(nx*ny,nt);  % Temperaturmatrix

f0 = zeros(nx, ny);
for i = 1 : nx
    for j = 1 : ny
        f0(i,j) = sin(x(i) / Lx * 2 * pi) * cos(y(j) / Ly * 2 * pi) + 1;
    end
end

[Y,X] = meshgrid(y,x);
mesh(X,Y,f0);
caxis([-0.5 2.5])
pbaspect([1 Ly/Lx 0.5])
setmesh('Anfangsbedingung','$x$','$y$','$f$','T_2D_inhomo_modal_inital_condition',0)

phi = zeros(N + 1, nx);  % Eigenfunktionen
psi = zeros(N + 1, ny);  % Eigenfunktionen
T = zeros((N + 1)^2, nt);  % Gewichtung

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

phi = phi';
psi = psi';

Psi = kron(psi,phi);

figure
for i = 1 : N + 1
    plot(x,phi(:,i))
    hold on
end
txt = ['$G = ',num2str(N),'$'];
TEXT = text(8,0.4,txt,'FontSize',60);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunktionen $\varphi$','$x$','','T_2D_inhomo_modal_Eigenfunctions_phi',0)

figure
for i = 1 : N + 1
    plot(y,psi(:,i))
    hold on
end
txt = ['$G = ',num2str(N),'$'];
TEXT = text(16,0.3,txt,'FontSize',60);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunktionen $\psi$','$y$','','T_2D_inhomo_modal_Eigenfunctions_psi',0)

u = zeros(nt, nx, ny);  % Anregung
Ua = zeros(N + 1, N + 1, nt);

u(:,round(3/dx+1),round(3/dy+1)) = 1 * sin(t - pi / 4);
u(:,round(5/dx+1),round(10/dy+1)) = -2 * sin(t);
u(:,round(7/dx+1),round(13/dy+1)) = 0.02 * t;
u(:,round(8/dx+1),round(18/dy+1)) = 0.05 * t;

for m = 1 : N + 1
    for n = 1 : N + 1
        for i = 1 : nt
            Ua(m,n,i) = u(i,round(3/dx+1),round(3/dy+1)) * phi(round(3/dx+1),m) * psi(round(3/dy+1),n) +...
                       u(i,round(5/dx+1),round(10/dy+1)) * phi(round(5/dx+1),m) * psi(round(10/dy+1),n) +...
                      u(i,round(7/dx+1),round(13/dy+1)) * phi(round(7/dx+1),m) * psi(round(13/dy+1),n) + ...
                       u(i,round(8/dx+1),round(18/dy+1)) * phi(round(8/dx+1),m) * psi(round(18/dy+1),n);
        end
    end
end

for i = 1 : nt
    U(:,i) = reshape(squeeze(Ua(:,:,i)),(N+1)^2,1);
end

T0 = zeros(N+1,N+1);
for i = 1 : N + 1
    for j = 1 : N + 1
        T0(i, j) = 0;
        for n = 1 : nx
            for m = 1 : ny
                T0(i, j) = T0(i, j) + f0(n,m) * phi(n,i) * psi(m,j) * dx * dy;
            end
        end
    end
end


kappa = zeros(N+1, N+1);
for i = 1 : N+1
    for j = 1 : N+1
        kappa(i,j) = 1 - dt * k * (lambda(i)^2 + sigma(j)^2);
    end
end

A = diag(reshape(kappa,(N+1)^2,1));

T(:,1) = reshape(T0,(N+1)^2,1);
for t = 2 : nt
    T(:,t) = A * T(:,t-1) + dt * U(:,t);
end

for t = 1 : nt
    f(:,t) = Psi * T(:,t);
    F(t,:,:) = reshape(f(:,t),nx,ny);

end


t = 0 : dt : Time;
figure
for n = 1 : 5 : nt
    mesh(X,Y,squeeze(F(n,:,:)))
    zlim([-1 3])
    caxis([-1 3])
    pbaspect([1 Ly/Lx 0.5])
    setmesh('','$x$','$y$','$f$','T_2D_inhomo_modal_2',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(n)),'$'];
    TEXT = text(8,0,0.5,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$G = ',num2str(N),'$'];
    TEXT = text(8,0,1.8,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    zticks([-1 1 3])
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    zticklabels({'-1','1','3'})
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
    mesh(X,Y,squeeze(F(a,:,:)))
    zlim([-1 3])
    caxis([-1 3])
    pbaspect([1 Ly/Lx 0.5])
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(a)),'$'];
    TEXT = text(8,0,0.5,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$G = ',num2str(N),'$'];
    TEXT = text(8,0,1.8,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    zticks([-1 1 3])
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    zticklabels({'-1','1','3'})
    drawnow
    name = ['T_2D_inhomo_modal_2_shot_',num2str(n)];
    setmesh('','$x$','$y$','$f$',name,printfigure)
end


close all

