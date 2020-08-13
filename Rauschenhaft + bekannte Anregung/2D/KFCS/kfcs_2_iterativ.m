% scr + kalman integrated + iterative

clc
clear
close all

printfigure = 0;

load('Messwerte_rh')
load basis_ml

ns = size(m_rh,1);

dctx = idct(eye(8));
dcty = idct(eye(8));
THETA = kron(dcty,dctx);
MT = THETA^-1;

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
[Y,X] = meshgrid(y,x);

t = 0 : Dt : 20;
nt = length(t);


% random select of sensors
number = 49;
S = zeros(nt,number);
S(:,1) = ceil(rand(nt,1)*ns);
for t = 1 : nt
    k = 2;
    while k <= number
        temp = ceil(rand()*ns);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end

% CS Basis
Phi_cs = zeros(nt, number, ns);
for t = 1 : nt
    for i = 1 : number
        Phi_cs(t,i,S(t,i)) = 1;
    end
end

z = zeros(ns, nt);
z(:,1) = MT * m_rh(:,1);


k = 0.1;


N = 7;  % Grad
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

index = p_index;
Phi_KF = zeros(ns, nx*ny);
for i = 1 : ns
    Phi_KF(i,index(i)) = 1;
end

H = Phi_KF * Psi;

T = zeros((N+1)^2,nt);
T(:,1) = (H' * H)^-1 * H' * m_rh(:,1);
Ce = eye((N+1)^2) * 1;
Cv = speye(ns) * 1;  % Messunsicherheit
Cw = speye((N+1)^2) * 0.1;  % Systemrauschen

f_e_integrated_2 = zeros(nx*ny,nt);
f_e_integrated_2(:,1) = Psi * T(:,1);

for t = 2 : nt
    h = squeeze(Phi_cs(t,:,:)) * THETA;
    y = squeeze(Phi_cs(t,:,:)) * m_rh();
    
    %dynamic weight
    Cv = eye(ns) * 5; 
    for i = 1 : size(S,2)
        Cv(S(t,i),S(t,i)) = 0.3;
    end
    
    Tp = A * T(:,t-1) ;
    Cp = A * Ce * A' + Cw;
    K = Cp * H' / (H * Cp * H' + Cv);


    % iterative
    zi = MT * Phi_KF * Psi * Tp;
    for iterative = 1 : 10
        e = 0.3;
        cvx_begin quiet;
            variable a(ns,1);
            minimize(norm(a,1));
            subject to;
                norm(h * (zi + a) - y) <= e;
        cvx_end;

        zi_new = zi + a;
        ypi = THETA * zi_new;
        
        Ti = (eye(size(K,1)) - K * H) * Tp + K * ypi;
        zi_kf = MT * Phi_KF * Psi * Ti;
        
        if norm(zi_kf - zi) < 100000000
            break
        else
            zi = zi_kf;
        end
    end
    z(:, t) = zi_kf;
    f_e_integrated_2(:,t) = Psi * Ti;

    ['t = ',num2str(t),',    iterativ = ',num2str(iterative)]
    Ce = (speye(size(K,1)) - K * H) * Cp;
end

F_e = zeros(nt,nx,ny);
for i = 1 : nt
    F_e(i,:,:) = reshape(f_e_integrated_2(:,i),nx,ny);
end

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
    txt = ['$M = ',num2str(number),'$'];
    TEXT = text(8,0,1.5,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,ns);
     name = ['Kalman_Modal_2D_',num2str(number),'_rh.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

f_rh_scr_kf_fdm_integrated_iterativ_2 = f_e_integrated_2;
save('f_rh_scr_kf_fdm_integrated_iterativ_2.mat','f_rh_scr_kf_fdm_integrated_iterativ_2')

