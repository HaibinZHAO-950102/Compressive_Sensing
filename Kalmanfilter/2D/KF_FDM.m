clc
clear
close all

printfigure = 1;

load('Messwerte_rh_2D')

M = 8;  % Anzahl der Messungen
Dt = 0.1; % time_step

dx = 2 * dx
x = 0 : dx : Lx;
nx = length(x);
dy = 2 * dy
y = 0 : dy : Ly;
ny = length(y);


Sx = round(linspace(1,8,M));  % benutzte Sensoren
Sy = round(linspace(1,8,M));  % benutzte Sensoren
m_rh = m_rh(1:Dt/dt:end,Sx,Sy);

px_index = round(linspace(1,nx,8));
py_index = round(linspace(1,ny,8))
;

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

Knoten = @(m,n) (n-1) * nx + m;
Punkte = @(knoten) [mod(knoten,nx) (knoten-mod(knoten,nx))/nx+1];

dt = 0.1;
A = k * dt / dx^2;
B = k * dt / dy^2;
Q = 1 + 2*A + 2*B;

listQi = [];
listnAi = [];
listnBi = [];
listQAi = [];
listQBi = [];
listQABi = [];

listQj = [];
listnAj = [];
listnBj = [];
listQAj = [];
listQBj = [];
listQABj = [];


for m = 1 : nx
    for n = 1 : ny
        listQi = [listQi Knoten(m,n)];
        listnAi = [listnAi Knoten(m,n)];
        listnAi = [listnAi Knoten(m,n)];
        listnBi = [listnBi Knoten(m,n)];
        listnBi = [listnBi Knoten(m,n)];
        
        listQj = [listQj Knoten(m,n)];
        listnAj = [listnAj Knoten(min(m+1,nx),n)];
        listnAj = [listnAj Knoten(max(m-1,1),n)];
        listnBj = [listnBj Knoten(m,min(n+1,ny))];
        listnBj = [listnBj Knoten(m,max(n-1,1))];
    end
end


invD =   sparse(listQi,listQj,Q,nx*ny,nx*ny)...
    + sparse(listnAi,listnAj,-A,nx*ny,nx*ny)...
    + sparse(listnBi,listnBj,-B,nx*ny,nx*ny);

% D = inv(invD);
   
load Anregung


Phi = zeros(M^2,nx*ny);
for i = 1 : M^2
    Phi(i,index(i)) = 1;
end
H = Phi;


phi = zeros(11, nx);  % Eigenfunktionen
psi = zeros(11, ny);  % Eigenfunktionen
T = zeros(11^2, nt);  % Gewichtung
N = 10;  % Grad
lambda = 0 : pi / Lx : N * pi / Lx;
sigma = 0 : pi / Ly : N * pi / Ly;

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


findexx = round(linspace(1,128,nx));
findexy = round(linspace(1,128,ny));
F0 = squeeze(F_sr(1,findexx,findexy));

f = zeros(nx*ny,nt);
f(:,1) = reshape(F0,nx*ny,1);

% for l = 2 : 201
%     f(:,l) = invD \ (f(:,l-1) + U(:,l-1));
% end

Ce = sparse(nx*ny, nx*ny);
Cv = speye(M^2) * 0.01;  % Messunsicherheit
Cw = Psi * eye(121)* Psi' * 5;

for t = 2 : nt
    t
    fp = invD \ f(:,t-1);
    Cp = invD \ Ce / invD' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);
    f(:,t) = (eye(size(K,1)) - K * H) * fp + K * reshape(m_rh(t,:,:),M^2,1);
    Ce = (eye(size(K,1)) - K * H) * Cp;
end

F_e = zeros(nt,nx,ny);
for t = 1 : nt
    F_e(t,:,:) = reshape(f(:,t),nx,ny);
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
     name = ['Kalman_FDM_2D_',num2str(M^2),'_rh.gif'];
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
    name = ['T_2D_FDM_64_shot_',num2str(n)];
    setmesh('','$x$','$y$','$f$',name,printfigure)
end

close all














