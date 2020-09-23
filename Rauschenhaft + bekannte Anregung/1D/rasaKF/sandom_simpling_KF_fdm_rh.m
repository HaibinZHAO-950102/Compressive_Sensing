clc
clear
close all

printfigure = 1;

load('Messwerte_rh')


Dt = 0.1;
m_rh = m_rh(:,1:Dt/dt:end);

N_time = size(m_rh,2);

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
D = speye(nx) - P*A;

u = zeros(nt, nx);
u(:,round(3/dx+1)) = 0.1 * sin(t - pi / 4) / dx;
u(:,round(5/dx+1)) = -0.2 * sin(t) / dx;
u(:,round(7/dx+1)) = 0.01 * t / dx;


Phi = zeros(N_time, number, nx);
for t = 1 : N_time
    for i = 1 : number
        Phi(t,i,p_index(S(t,i))) = 1;
    end
end

f_e = zeros(nx, nt);
f_e(:,1) = f_mu(:,1);
Ce = speye(nx) * 1;
Cv = speye(nx) * 0.1;  % Messunsicherheit
Cw = speye(nx) * 0.1;  % Systemrauschen


for t = 2 : N_time
    h = squeeze(Phi(t,:,:));
    cv = squeeze(Phi(t,:,:)) * Cv * squeeze(Phi(t,:,:))';
    y = m_rh(S(t,:),t);
    fp = D \ f_e(:, t-1) + D^-1 * u(t-1,:)' * Dt;
    Cp = D \ (Ce + Cw) / D';
    K = Cp * h' / (h * Cp * h' + cv);
    f_e(:,t) = (eye(size(K,1)) - K * h) * fp + K * y;
    Ce = (eye(size(K,1)) - K * h) * Cp;
end


figure
for n = 1 : 0.5/Dt : N_time
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    clf
    plot(x, f_sr(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e(:,n),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(n,:))),m_rh(S(n,:),n),'r.','Markersize',40)
    legend('Signal with Noise','Estimated Signal','Measure Points')
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
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['random_simpling_KF_fdm_rh.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

f_rh_random_simpling_KF_fdm = f_e;
save('f_rh_random_simpling_KF_fdm.mat','f_rh_random_simpling_KF_fdm')


