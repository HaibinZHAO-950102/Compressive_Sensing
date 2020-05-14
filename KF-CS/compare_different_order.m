clc
clear
close all
printfigure = 0;

load('FDM')
load('Modal_10')
load('Modal_50')
load('Modal_200')
load('Modal_300')

Time = 20;
t = 0 : 0.1 : Time;
x = 0 : 0.1 : 10;

for n = 1 : 201
    clf
    plot(x, f_fdm(n,:),'k-','LineWidth',1)
    hold on
    plot(x, f_modal_10(n,:),'r-','LineWidth',1)
    hold on
    plot(x, f_modal_50(n,:),'g-','LineWidth',1)
    hold on
    plot(x, f_modal_200(n,:),'b-','LineWidth',1)
    hold on
    plot(x, f_modal_300(n,:),'m-','LineWidth',1)
    hold on
    legend('FDM','Modal N=10','Modal N=50','Modal N=200','Modal N=300')
    xlim([0 10])
    ylim([0 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(Time/200*(n-1)),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'Compare_N.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'Compare_N.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

figure
plot(x, f_fdm(end,:),'k-','LineWidth',1)
hold on
plot(x, f_modal_10(end,:),'r-','LineWidth',1)
hold on
plot(x, f_modal_50(end,:),'g-','LineWidth',1)
hold on
plot(x, f_modal_200(end,:),'b-','LineWidth',1)
hold on
plot(x, f_modal_300(end,:),'m-','LineWidth',1)
hold on
legend('FDM','Modal N=10','Modal N=50','Modal N=200','Modal N=300')
xlim([0 10])
ylim([0 2.5])
set(gca,'Fontsize',20)
set(gca,'fontname','times new Roman')
T = title('Temperature Distribution','fontsize',40);
set(T,'Interpreter','latex')
T = xlabel('$x$','fontsize',30);
set(T,'Interpreter','latex')
T = ylabel('$T$','fontsize',30);
set(T,'Interpreter','latex')
set(gcf,'outerposition',get(0,'screensize'));
txt = ['$t = 20$'];
T = text(0.8,0.6,txt,'FontSize',30);
set(T,'Interpreter','latex')
drawnow
setplt('Compare different order','$x$','$T$','Compare_different_order',printfigure)

