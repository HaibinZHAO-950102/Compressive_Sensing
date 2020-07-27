clc
clear
close all

printfigure = 1;

load('Messwerte')

M = 36;  % Anzahl der Messungen
Dt = 0.1; % time_step

S = round(linspace(1,64,M));  % benutzte Sensoren
m = m(S,1:Dt/dt:end);


x = 0 : dx : 10;
nx = length(x);
t = 0 : Dt : 20;
nt = length(t);

k = 0.1;
E = sparse(2:nx,1:nx-1,1,nx,nx);
P = Dt*k/dx^2;
A = -2*speye(nx)+(E+E'); 
A(1,1) = -1; % Neumann Boundary
A(nx,nx) = -1; % Neumann Boundary
U = speye(nx) - P*A;


N_time = size(m,2);


Phi = zeros(length(S),nx);
for i = 1 : length(S)
    Phi(i,p_index(S(i))) = 1;
end

H = Phi;


f_e = zeros(nx, nt);
Ce = speye(nx) * 1;
Cv = speye(length(S)) * 0.1;  % Messunsicherheit
Cw = speye(nx) * 1;  % Systemrauschen
for t = 2 : N_time
    fp = U \ f_e(:, t-1);
    Cp = U \ (Ce + Cw) / U';
    K = Cp * H' / (H * Cp * H' + Cv);
    f_e(:,t) = (eye(size(K,1)) - K * H) * fp + K * m(:,t);
    Ce = (speye(size(K,1)) - K * H) * Cp;
end


figure
for n = 1 : 0.5/Dt : N_time
    clf
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    plot(x, f(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e(:,n),'c-','LineWidth',5)
    hold on
    plot(p(S),m(:,n),'r.','Markersize',40)
    legend('Signal','Signal Estimated','Measure Points')
    xlim([0 10])
    ylim([-0.5 2.5])
    setplt('Temperature Distribution','$x$','$f$','Temperature Distribution',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(M),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['Kalman_fdm_',num2str(M),'.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% f_e_kalman_fdm_1D_40 = f_e;
% save('f_e_kalman_fdm_1D_40.mat','f_e_kalman_fdm_1D_40')

