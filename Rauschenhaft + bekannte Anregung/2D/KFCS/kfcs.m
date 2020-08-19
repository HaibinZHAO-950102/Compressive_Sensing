clc
clear
close all

printfigure = 1;

load('Messwerte_rh')
load CS_SCR_L1

m_rh = yp;

ns = size(m_rh,1);
    
M = 8;  % Anzahl der Messungen
Dt = 0.1; % time_step

dt = 0.1;
Dt = 0.1; % time_step

dx = 10/127;
dy = 20/127;
Lx = 10;
Ly = 20;

x = 0 : dx : Lx;
nx = length(x);
y = 0 : dy : Ly;
ny = length(y);


Sensor = zeros(nx*ny,1);
Sensor(p_index) = 1;
index = p_index;

x = 0 : dx : 10;
nx = length(x);
y = 0 : dy : 20;
ny = length(y);
t = 0 : Dt : 20;
nt = length(t);

k = 0.1;  % Waermeleitfaehigkeit in cm^2/s

N = min(10,M-1);  % Grad
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


Phi = zeros(ns,nx*ny);
for i = 1 : ns
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
    U(:,i) = reshape(squeeze(Ua(:,:,i)),ns,1);
end


T = zeros((N+1)^2,nt);
T(:,1) = (Psi' * Psi)^-1 * Psi' * f(:,1);
Ce = zeros((N+1)^2, (N+1)^2, nt);
Ce(:,:,1) = eye((N+1)^2) * 1;
Cv = speye(ns) * 2;  % Messunsicherheit
Cw = speye((N+1)^2) * 1e-4;  % Systemrauschen

for t = 2 : nt
    Tp = A * T(:,t-1) + Dt * U(:,t-1);
    Cp = A * Ce(:,:, t-1) * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    T(:,t) = (eye(size(K,1)) - K * H) * Tp + K * m_rh(:,t);
    Ce(:,:,t) = (eye(size(K,1)) - K * H) * Cp;
end

f_e = zeros(nx*ny,t);
F_e = zeros(nt,nx,ny);
for t = 1 : nt
    f_e(:,t) = Psi * T(:,t);
    F_e(t,:,:) = reshape(f_e(:,t),nx,ny);
end

index = 1:nx*ny;
figure
for n = 1 : 0.5/Dt : nt
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    clf
    plot(index, f(:,t_f),'k-','LineWidth',1)
    hold on
    plot(index, f_e(:,n),'b-','LineWidth',1)
    hold on
    plot(index(p_index(S(n,:))),m_rh(S(n,:),n)','r.','Markersize',20)
    legend('Signal','Signal Estimated','Measure Points')
    xlim([0 nx*ny])
    ylim([-0.5 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    z = title('Temperature Distribution','fontsize',40);
    set(z,'Interpreter','latex')
    z = xlabel('$x$','fontsize',30);
    set(z,'Interpreter','latex')
    z = ylabel('$f$','fontsize',30);
    set(z,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    z = text(0.8,0.6,txt,'FontSize',30);
    set(z,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['KFCS_1D.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end




[Y,X] = meshgrid(y,x);
t = 0 : Dt : 20;
figure
for n = 1 : 5 : nt
    mesh(X,Y,squeeze(F_e(n,:,:)))
    zlim([-1 3])
    caxis([-1 3])
    pbaspect([1 Ly/Lx 0.5])
    setmesh('Temperature Distribution','$x$','$y$','$f$','T_2D_inhomo_modal_2',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(n)),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    txt = ['$M = 25$'];
    TEXT = text(8,0,1.5,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
     name = ['KFCS_2D.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

f_e_kfcs = f_e;
save('f_e_kfcs.mat','f_e_kfcs')
