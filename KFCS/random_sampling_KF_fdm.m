clc
clear
close all

printfigure = 1;

load('Messwerte')


Dt = 0.1;
m = m(:,1:Dt/dt:end);

N_time = size(m,2);

number = 12;
S = zeros(N_time,number);
S(:,1) = ceil(rand(N_time,1)*64);
for t = 1 : N_time
    k = 2;
    while k <= number
        temp = ceil(rand()*64);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end

M = number;

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



Phi = zeros(N_time, number, nx);
for t = 1 : N_time
    for i = 1 : number
        Phi(t,i,p_index(S(t,i))) = 1;
    end
end

f_e = zeros(nx, nt);
f_e(:,1) = f(:,1);
Ce = speye(nx) * 1;
Cv = speye(nx) * 0.1;  % Messunsicherheit
Cw = speye(nx) * 10;  % Systemrauschen


for t = 2 : N_time
    h = squeeze(Phi(t,:,:));
    cv = squeeze(Phi(t,:,:)) * Cv * squeeze(Phi(t,:,:))';
    y = m(S(t,:),t);
    fp = U \ f_e(:, t-1);
    Cp = U \ (Ce + Cw) / U';
    K = Cp * h' / (h * Cp * h' + cv);
    f_e(:,t) = (eye(size(K,1)) - K * h) * fp + K * y;
    Ce = (eye(size(K,1)) - K * h) * Cp;
end


figure
for n = 1 : 0.5/Dt : N_time
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    clf
    plot(x, f(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e(:,n),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(n,:))),m(S(n,:),n),'r.','Markersize',40)
    legend('Signal','geschätztes Signal','Messungen')
    xlim([0 10])
    ylim([-0.5 2.5])
    setplt('Temperaturverteilung','$x$','$f$','Temperaturverteilung',0)
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(0.8,0.4,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['random_simpling_KF_fdm.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

t = 0 : 0.1:20;
for n = 0 : 5
    figure
    timef = 2000 / 5 * n + 1
    timefe = 200 / 5 * n + 1
    plot(x, f(:,timef),'k-','LineWidth',5)
    hold on
    plot(x, f_e(:,timefe),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(timefe,:))),m(S(timefe,:),timefe),'r.','Markersize',40)
    legend('Signal','geschätztes Signal','Messungen')
    xlim([0 10])
    ylim([-0.5 2.5])
    txt = ['$t = ',num2str(4*n),'$'];
    T = text(0.8,0.4,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['RSKF_FDM_shot_',num2str(n+1)];
    setplt('','$x$','$f$',name,printfigure)
end

% f_random_sampling_KF_fdm = f_e;
% save('f_random_sampling_KF_fdm.mat','f_random_sampling_KF_fdm')

close all
