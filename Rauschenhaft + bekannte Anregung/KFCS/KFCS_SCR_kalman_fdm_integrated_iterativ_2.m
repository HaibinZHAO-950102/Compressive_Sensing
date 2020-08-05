% mt + kalman integrated + iterative

clc
clear
close all

printfigure = 0;

load('Messwerte_rh')
load basis_ml

x = 0 : dx : 10;
nx = length(x);

Dt = 0.1;
m_rh = m_rh(:,1:Dt/dt:end);

nt = size(m_rh,2);

% random select of sensors
number = 12;
S = zeros(nt,number);
S(:,1) = ceil(rand(nt,1)*64);
for t = 1 : nt
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


% CS Basis
Phi_cs = zeros(nt, number, 64);
for t = 1 : nt
    for i = 1 : number
        Phi_cs(t,i,S(t,i)) = 1;
    end
end

z = zeros(64, nt);
z(:,1) = MT * m_rh(:,1);


k = 0.1;


E = sparse(2:nx,1:nx-1,1,nx,nx);
P = Dt*k/dx^2;
A = -2*speye(nx)+(E+E'); 
A(1,1) = -1; % Neumann Boundary
A(nx,nx) = -1; % Neumann Boundary
D = speye(nx) - P*A;

t = 0 : Dt : 20;
u = zeros(nt, nx);
u(:,round(3/dx+1)) = 0.1 * sin(t - pi / 4) / dx;
u(:,round(5/dx+1)) = -0.2 * sin(t) / dx;
u(:,round(7/dx+1)) = 0.01 * t / dx;


Phi_KF = zeros(length(p_index),nx);
for i = 1 : length(p_index)
    Phi_KF(i,p_index(i)) = 1;
end


H = zeros(64, nx);
for i = 1 : 64
    H(i,p_index(i)) = 1;
end



f_e_integrated_2 = zeros(nx, nt);
f_e_integrated_2(:,1) = f_mu(:,1);
Ce = speye(nx) * 1;
Cv = speye(length(p_index)) * 1;  % Messunsicherheit
Cw = speye(nx) * 0.1;  % Systemrauschen
for t = 2 : nt
    h = squeeze(Phi_cs(t,:,:)) * THETA;
    y = m_rh(S(t,:),t);
    
    %dynamic weight
    Cv = eye(length(p_index)) * 1; 
    for i = 1 : size(S,2)
        Cv(S(t,i),S(t,i)) = 0.01;
    end
    
    fp = D \ f_e_integrated_2(:, t-1) + D^-1 * u(t-1,:)' * Dt;
    Cp = D \ (Ce + Cw) / D';
    K = Cp * H' / (H * Cp * H' + Cv);


    % iterative
    zi = MT * Phi_KF * fp;
    for iterative = 1 : 10
        e = 0.01;
        cvx_begin quiet;
            variable a(64,1);
            minimize(norm(a,1));
            subject to;
                norm(h * (zi + a) - y) <= e;
        cvx_end;

        zi_new = zi + a;
        ypi = THETA * zi_new;
        
        fi = (eye(size(K,1)) - K * H) * fp + K * ypi;
        zi_kf = MT * Phi_KF * fi;
        
        if norm(zi_kf - zi) < 10
            break
        else
            zi = zi_kf;
        end
    end
    z(:, t) = zi_kf;
    f_e_integrated_2(:,t) = fi;

    ['t = ',num2str(t),',    iterativ = ',num2str(iterative)]
    Ce = (speye(size(K,1)) - K * H) * Cp;
end

figure
for n = 1 : 0.5/Dt : nt
    clf
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    plot(x, f_sr(:,t_f),'k-','LineWidth',5)
    hold on
    plot(x, f_e_integrated_2(:,n),'c-','LineWidth',5)
    hold on
    plot(x(p_index(S(n,:))),m_rh(S(n,:),n),'r.','Markersize',40)
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
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_SCR_KF_FDM_integrated_iterativ_rh_2.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% f_rh_scr_kf_fdm_integrated_iterativ_2 = f_e_integrated_2;
% save('f_rh_scr_kf_fdm_integrated_iterativ_2.mat','f_rh_scr_kf_fdm_integrated_iterativ_2')

