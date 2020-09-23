% mt cs
clc
clear
close all

printfigure = 1;

load('CS_SCR_L1');
load('Messwerte');

x = 0 : dx : 10;
nx = length(x);
Dt = 0.1;
t = 0 : Dt : 20;
nt = length(t);

f = f(:,1:Dt/dt:end);

m = yp;
xm = x(p_index);

figure
for n = 1 : 0.5/Dt : nt
    clf
    Signal = plot(x,f(:,n),'k-','linewidth',5)
    set(Signal,'Color',[0.9 0.9 0.9])
    hold on
    PM = plot(xm,m(:,n),'.','Markersize',30);
    set(PM,'Color',[0.1961 0.6314 0.5373])
    hold on
    MS = plot(x(p_index(S(n,:))),m(S(n,:),n),'r.','Markersize',40)
    legend('Signal','Pseudo Measure Points','Measurements')
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
    txt = ['$M = 12$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['Pseudo Measurements SCR.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end
