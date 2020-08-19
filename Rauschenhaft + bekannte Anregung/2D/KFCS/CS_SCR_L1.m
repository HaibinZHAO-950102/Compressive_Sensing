% meine approximation

clc
clear
close all

printfigure = 0;

load('Messwerte_rh')
load OSB
ns = size(m_rh,1);

% dctx = idct(eye(8));
% dcty = idct(eye(8));
% THETA = kron(dcty,dctx);
% ST = THETA^-1;

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
number = 25;
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
Phi = zeros(nt, number, ns);
for t = 1 : nt
    for i = 1 : number
        Phi(t,i,S(t,i)) = 1;
    end
end

z = zeros(ns, nt);
z(:,1) = ST * f(p_index,1);

for t = 2 : nt
    t
    h = squeeze(Phi(t,:,:)) * THETA;
    y = m_rh(S(t,:),t);
    
    e = 0.01;
    cvx_begin quiet;
        variable a(ns,1);
        minimize(norm(a,1));
        subject to;
            norm(h * (z(:,t-1) + a) - y) <= e;
    cvx_end;
    z(:, t) = z(:,t-1) + a;
end

yp = zeros(ns,nt);
for n = 1 : nt
    yp(:,n) = THETA * z(:,n);
end

index = 1 : nx * ny;

figure
for n = 1 : 0.5/Dt : nt
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    clf
%     plot(index, f(:,t_f),'k-','LineWidth',2)
%     hold on
    plot(index(p_index), m_rh(:,n),'g.','markersize',20)
    hold on
    plot(index(p_index), yp(:,n),'b.','markersize',20)
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
    txt = ['$M = ',num2str(number),'$'];
    z = text(0.8,0.8,txt,'FontSize',30);
    set(z,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['CS_SCR_L1.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

save('CS_SCR_L1.mat','yp','S')




