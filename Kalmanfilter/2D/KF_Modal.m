clc
clear
close all

printfigure = 1;

load('Messwerte_rh_2D')

M = 8;  % Anzahl der Messungen
Dt = 0.1; % time_step

x = 0 : dx : Lx;
nx = length(x);
y = 0 : dy : Ly;
ny = length(y);

Sx = round(linspace(1,8,M));  % benutzte Sensoren
Sy = round(linspace(1,8,M));  % benutzte Sensoren
m_rh = m_rh(1:Dt/dt:end,Sx,Sy);

Sensor = zeros(nx,ny);
for i = 1 : M
    for j = 1 : M
        Sensor(px_index(Sx(i)),py_index(Sy(j))) = 1;
    end
end
Sensor = reshape(Sensor,nx*ny,1);
index = find(Sensor == 1);

x = 0 : dx : 10;
nx = length(x);
y = 0 : dy : 20;
ny = length(y);
t = 0 : Dt : 20;
nt = length(t);

k = 0.1;  % Waermeleitfaehigkeit in cm^2/s

N = M-1;  % Grad
lambda = 0 : pi / Lx : N * pi / Lx;
sigma = 0 : pi / Ly : N * pi / Ly;


kappa = zeros(N+1, N+1);
for i = 1 : N+1
    for j = 1 : N+1
        kappa(i,j) = 1 - dt * k * (lambda(i)^2 + sigma(j)^2);
    end
end

A = diag(reshape(kappa,(N+1)^2,1));



phi = zeros(N + 1, nx);  % Eigenfunktionen
psi = zeros(N + 1, ny);  % Eigenfunktionen

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


Phi = zeros(M^2,nx*ny);
for i = 1 : M^2
    Phi(i,index(i)) = 1;
end
H = Phi * Psi;

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


T = zeros((N+1)^2,nt);
T(:,1) = (Psi' * Psi)^-1 * Psi' * reshape(F_sr(1,:,:),128^2,1);
Ce = zeros((N+1)^2, (N+1)^2, nt);
Ce(:,:,1) = eye((N+1)^2) * 1;
Cv = speye(M^2) * 0.1;  % Messunsicherheit
Cw = speye((N+1)^2) * 1e-3;  % Systemrauschen

for t = 2 : nt
    Tp = A * T(:,t-1) + Dt * U(:,t-1);
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,t) = (eye(size(K,1)) - K * H) * Tp + K * reshape(m_rh(t,:,:),M^2,1);
    Ce(:,:,t) = (eye(size(K,1)) - K * H) * Cp;
end

f_e = zeros(nx*ny,t);
F_e = zeros(nt,nx,ny);
for t = 1 : nt
    f_e(:,t) = Psi * T(:,t);
    F_e(t,:,:) = reshape(f_e(:,t),nx,ny);
end


[Y,X] = meshgrid(y,x);
t = 0 : Dt : 20;
figure
for n = 1 : 5 : nt
    mesh(X,Y,squeeze(F_e(n,:,:)))
    zlim([-1 3])
    caxis([-1 3])
    pbaspect([1 Ly/Lx 0.5])
    setmesh('','$x$','$y$','$f$','T_2D_inhomo_modal_2',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(n)),'$'];
    TEXT = text(8,0,0.5,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$M = ',num2str(M^2),'$'];
    TEXT = text(8,0,1.8,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    zticks([-1 1 3])
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    zticks([-1 1 3])
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
     name = ['Kalman_Modal_2D_',num2str(M^2),'_rh.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% f_e_kf = f_e;
% save('f_e_kf.mat','f_e_kf')

t = 0 : 0.1:20;
[Y,X] = meshgrid(y,x);

for n = 0 : 5
    timefe = 200 / 5 * n + 1;

    figure
    mesh(X,Y,squeeze(F_e(timefe,:,:)))
    zlim([-1 3])
    caxis([-1 3])
    pbaspect([1 Ly/Lx 0.5])
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(timefe)),'$'];
    TEXT = text(8,0,0.5,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    txt = ['$M = ',num2str(64),'$'];
    TEXT = text(8,0,1.8,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    zticks([-1 1 3])
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    zticks([-1 1 3])
    drawnow
    name = ['T_2D_modal_64_shot_',num2str(n)];
    setmesh('','$x$','$y$','$f$',name,printfigure)
end

close all


