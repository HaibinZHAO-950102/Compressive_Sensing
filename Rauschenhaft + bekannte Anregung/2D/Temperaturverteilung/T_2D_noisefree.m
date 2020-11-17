clc
clear
close all

printfigure = 0;

Lx = 10;  
Ly = 20;  
Time = 20;
dx = Lx / 63;
dy = Ly / 63;
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
setmesh('Initial Condition','$x$','$y$','$f$','T_2D_inhomo_modal_inital_condition',0)

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
TEXT = text(8,0.4,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunctions $\varphi$','$x$','$value$','T_2D_inhomo_modal_Eigenfunctions_phi',0)

figure
for i = 1 : N + 1
    plot(y,psi(:,i))
    hold on
end
txt = ['$G = ',num2str(N),'$'];
TEXT = text(16,0.3,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunctions $\psi$','$y$','$value$','T_2D_inhomo_modal_Eigenfunctions_psi',0)

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

w = randn((N+1)^2, nt) * 0.02;
v = randn(nx*ny, nt) * 0.05;
T(:,1) = reshape(T0,(N+1)^2,1);
for t = 2 : nt
    T(:,t) = A * T(:,t-1) + dt * U(:,t);
end

for t = 1 : nt
    f(:,t) = Psi * T(:,t);
    V(t,:,:) = reshape(v(:,t),nx,ny);
    F(t,:,:) = reshape(f(:,t),nx,ny);

end

F_mu = F;

t = 0 : dt : Time;
figure
for n = 1 : 5 : nt
    mesh(X,Y,squeeze(F(n,:,:)))
    zlim([-1 3])
    caxis([-1 3])
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
             imwrite(imind,cm,'T_2D_sr.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'T_2D_sr.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end


figure
for n = 1 : 5 : nt
    mesh(X,Y,squeeze(F_mu(n,:,:)))
    zlim([-1 3])
    caxis([-1 3])
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
             imwrite(imind,cm,'T_2D_mu.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'T_2D_mu.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



% figure
% for n = 1 : 6
%     a = round((nt - 1) / 5 * (n - 1)) + 1;
%     mesh(X,Y,squeeze(F(a,:,:)))
%     zlim([-1 3])
%     caxis([-1 3])
%     pbaspect([1 Ly/Lx 0.5])
%     set(gcf,'outerposition',get(0,'screensize'));
%     txt = ['$t = ',num2str(t(a)),'$'];
%     TEXT = text(8,0,0.6,txt,'FontSize',30);
%     set(TEXT,'Interpreter','latex')
%     txt = ['$G = ',num2str(N),'$'];
%     TEXT = text(8,0,1.5,txt,'FontSize',30);
%     set(TEXT,'Interpreter','latex')
%     drawnow
%     name = ['T_2D_sr_shot_',num2str(n)];
%     setmesh('Temperature Distribution','$x$','$y$','$f$',name,printfigure)
% end
% 
% figure
% for n = 1 : 6
%     a = round((nt - 1) / 5 * (n - 1)) + 1;
%     mesh(X,Y,squeeze(F_mu(a,:,:)))
%     zlim([-1 3])
%     caxis([-1 3])
%     pbaspect([1 Ly/Lx 0.5])
%     set(gcf,'outerposition',get(0,'screensize'));
%     txt = ['$t = ',num2str(t(a)),'$'];
%     TEXT = text(8,0,0.6,txt,'FontSize',30);
%     set(TEXT,'Interpreter','latex')
%     txt = ['$G = ',num2str(N),'$'];
%     TEXT = text(8,0,1.5,txt,'FontSize',30);
%     set(TEXT,'Interpreter','latex')
%     drawnow
%     name = ['T_2D_mu_shot_',num2str(n)];
%     setmesh('Temperature Distribution','$x$','$y$','$f$',name,printfigure)
% end

% Measurement
% positions
px_index = round(linspace(1,nx,8)); 
py_index = round(linspace(1,ny,8));

px = x(px_index);
py = y(py_index);


% m_rh = F_mu(:,px_index,py_index);   
% F_sr = F;
save('T_2D_noisefree','F')

