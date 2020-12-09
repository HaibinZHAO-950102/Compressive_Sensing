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
    plot(x,f(:,n),'k-','linewidth',5)
    hold on
    PM = plot(xm,m(:,n),'.','Markersize',40);
    set(PM,'Color',[0.1961 0.6314 0.5373])
    hold on
    MS = plot(x(p_index(S(n,:))),m(S(n,:),n),'r.','Markersize',40)
    legend('Signal','Pseudo-Messungen','echte Messungen')
    xlim([0 10])
    ylim([-0.5 2.5])
    setplt('','$x$','$f$','',0);
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(0.8,0.4,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    txt = ['$M = 12$'];
    T = text(0.8,0.8,txt,'FontSize',60);
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

t = 0 : 0.1:20;
for n = 0 : 5
    figure
    timefe = 200 / 5 * n + 1
    plot(x,f(:,timefe),'k-','linewidth',5)
    hold on
    PM = plot(xm,m(:,timefe),'.','Markersize',40);
    set(PM,'Color',[0.1961 0.6314 0.5373])
    hold on
    MS = plot(x(p_index(S(timefe,:))),m(S(timefe,:),timefe),'r.','Markersize',40)
    legend('Signal','Pseudo-Messungen','echte Messungen')

    xlim([0 10])
    ylim([-0.5 2.5])
    txt = ['$t = ',num2str(4*n),'$'];
    T = text(0.8,0.4,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    txt = ['$M = 12$'];
    T = text(0.8,0.8,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['Pseudo_Messungen_shot_',num2str(n+1)];
    setplt('','$x$','$f$',name,printfigure)
end


close all
