% meine approximation

clc
clear
close all

printfigure = 1;

load('Messwerte')
load basis_ml

x = 0 : dx : 10;
nx = length(x);
Dt = 0.1;
m = m(:,1:Dt/dt:end);
nt = size(m,2);

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



Phi = zeros(nt, number, 64);
for t = 1 : nt
    for i = 1 : number
        Phi(t,i,S(t,i)) = 1;
    end
end

z = zeros(64, nt);
z(:,1) = MT * f(p_index,1);

% for t = 2 : nt
%     t
%     h = squeeze(Phi(t,:,:)) * THETA;
%     y = m(S(t,:),t);
%     
%     e = 0.03;
%     cvx_begin quiet;
%         variable a(64,1);
%         minimize(norm(a,1));
%         subject to;
%             norm(h * (z(:,t-1) + a) - y) <= e;
%     cvx_end;
%     z(:, t) = z(:,t-1) + a;
% end
% 
% yp = zeros(64,nt);
% for n = 1 : nt
%     yp(:,n) = THETA * z(:,n);
% end
% 

figure
for n = 1 : 0.5/Dt : nt
    t = (n-1) * Dt;
    t_f = t/dt + 1;
    clf
    Signal = plot(x, f(:,t_f),'-','LineWidth',5)
    set(Signal,'Color',[0.9 0.9 0.9])
    hold on
    MS = plot(x(p_index),m(:,n),'k.','Markersize',20)
    set(MS,'Color',[0.5 0.5 0.5])
    hold on
    plot(x(p_index(S(n,:))),m(S(n,:),n)','r.','Markersize',40)
    legend('Signal','possible Sensor Positions','real Measure Points')
    xlim([0 10])
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
    name = ['zufall_Messung.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% save('CS_SCR_L1.mat','yp','S')




