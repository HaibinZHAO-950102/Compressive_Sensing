clc
clear
close all

printfigure = 0;

load('f_fdm_1')
load('f_modal_1')
load('f_fdm_2')
load('f_modal_2')

Length = 10;  % Stablaenge
Time = 5;   % Zetiraum

figure
for n = 1 : 51
    t_fdm = 1 * (n - 1) + 1;
    t_modal = 10 * (n - 1) + 1;
    clf
    plot(x_fdm_1, f_fdm_1(t_fdm,:),'b-','LineWidth',5)
    hold on
    plot(x_modal_1, f_modal_1(t_modal,:),'r-','LineWidth',5)
    hold on
    legend('FDM','Modal')
    xlim([0 10])
    ylim([0 2])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(Time/50*(n-1)),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'TV_modal_fdm_1.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'TV_modal_fdm_1.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end


Length = 10;  % Stablaenge
Time = 20;   % Zetiraum

figure
for n = 1 : 51
    t_fdm = 4 * (n - 1) + 1;
    t_modal = 40 * (n - 1) + 1;
    clf
    plot(x_fdm_2, f_fdm_2(t_fdm,:),'b-','LineWidth',5)
    hold on
    plot(x_modal_2, f_modal_2(t_modal,:),'r-','LineWidth',5)
    hold on
    legend('FDM','Modal')
    xlim([0 10])
    ylim([0 2])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(Time/50*(n-1)),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'TV_modal_fdm_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'TV_modal_fdm_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end
