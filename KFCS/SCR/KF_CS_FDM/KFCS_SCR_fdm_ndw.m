% mt + kalman

clc
clear
close all

printfigure = 1;

load('CS_SCR_L1');
load('Messwerte')

x = 0 : dx : 10;
nx = length(x);
Dt = 0.1;
t = 0 : Dt : 20;
nt = length(t);


m = yp;
xm = x(p_index);

k = 0.1;
E = sparse(2:nx,1:nx-1,1,nx,nx);
P = Dt*k/dx^2;
A = -2*speye(nx)+(E+E'); 
A(1,1) = -1; % Neumann Boundary
A(nx,nx) = -1; % Neumann Boundary
U = speye(nx) - P*A;


N_time = size(m,2);


Phi = zeros(length(p_index),nx);
for i = 1 : length(p_index)
    Phi(i,p_index(i)) = 1;
end

H = Phi;

f_e = zeros(nx, nt);
f_e = f(:,1);
Ce = speye(nx) * 1;
Cv = speye(length(S)) * 1;  % Messunsicherheit
Cw = speye(nx) * 5;  % Systemrauschen
for n = 2 : nt
    n
    Cv = eye(length(p_index)) * 1; 
    
    fp = U \ f_e(:, n-1);
    Cp = U \ (Ce + Cw) / U';
    K = Cp * H' / (H * Cp * H' + Cv);
    f_e(:,n) = (eye(size(K,1)) - K * H) * fp + K * m(:,n);
    Ce = (speye(size(K,1)) - K * H) * Cp;
end


figure
for n = 1 : 0.5/Dt : nt
    clf
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    plot(x, f(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e(:,n),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(n,:))),m(S(n,:),n),'r.','Markersize',40)
    legend('Signal','Signal Estimated','Measure Poins')
    xlim([0 10])
    ylim([-0.5 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$f$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(size(S,2)),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_SCR_KF_fdm_ndw.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% f_e_scr_kf_fdm_ndw = f_e;
% save('f_e_scr_kf_fdm_ndw.mat','f_e_scr_kf_fdm_ndw')