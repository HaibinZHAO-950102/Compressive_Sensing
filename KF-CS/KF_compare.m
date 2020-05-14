clc
clear
close all
printfigure = 1;

load('KFCS')
load('f_e_kalman_6')
load('Messwerte')

x = 0 : 0.01:10;

figure
for n = 1 : 400 : size(f_e,1)
    clf
    plot(x, f(n,:),'k-','LineWidth',5)
    hold on
    plot(x, f_e(n,:),'r-','LineWidth',5)
    hold on
    plot(x, f_e_kalman_6(n,:),'b-','LineWidth',5)
    hold on
    legend('Origin','KF-CS','KF')
    xlim([0 10])
    ylim([-0.5 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*0.001),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = 6$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['compare KFCS and CS.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



error_kfcs = f - f_e;
error_cs = f - f_e_kalman_6;
zero = zeros(size(f));
figure
for n = 1 : 400 : size(f_e,1)
    clf
    plot(x, error_kfcs(n,:),'r-','LineWidth',5)
    hold on
    plot(x, error_cs(n,:),'b-','LineWidth',5)
    hold on
    hold on
    plot(x, zero(n,:),'k--','LineWidth',5)
    legend('KF-CS','KF')
    xlim([0 10])
    ylim([-1 1])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Error','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*0.001),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = 6$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['compare KFCS and CS error.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end