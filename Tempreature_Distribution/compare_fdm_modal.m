clc
clear
close all

printfigure = 0;

load('f_fdm_homo_1')
load('f_fdm_homo_2')
load('f_fdm_inhomo_1')
load('f_fdm_inhomo_2')
load('f_modal_homo_1')
load('f_modal_homo_2')
load('f_modal_inhomo_1')
load('f_modal_inhomo_2')

% homo 1

Length = 10;  % Stablaenge
Time = 200;   % Zetiraum

figure
for n = 1 : 51
    t_fdm = Time / 1 / 50 * (n - 1) + 1;
    t_modal = Time / 0.1 / 50 * (n - 1) + 1;
    clf
    plot(x_fdm_homo_1, f_fdm_homo_1(t_fdm,:),'b-','LineWidth',5)
    hold on
    plot(x_modal_homo_1, f_modal_homo_1(t_modal,:),'r-','LineWidth',5)
    hold on
    legend('FDM','Modal')
    xlim([0 10])
    ylim([-1 1])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$f$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(Time/50*(n-1)),'$'];
    T = text(8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'Compare_MA_FDM_homo_1.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'Compare_MA_FDM_homo_1.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

% homo 2

Length = 10;  % Stablaenge
Time = 0.5;   % Zetiraum

figure
for n = 1 : 51
    t_fdm = Time / 50 / 0.01 * (n - 1) + 1;
    t_modal = Time / 50 / 5e-4 * (n - 1) + 1;
    clf
    plot(x_fdm_homo_2, f_fdm_homo_2(t_fdm,:),'b-','LineWidth',5)
    hold on
    plot(x_modal_homo_2, f_modal_homo_2(t_modal,:),'r-','LineWidth',5)
    hold on
    legend('FDM','Modal')
    xlim([0 10])
    ylim([0 1])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$f$','fontsize',30);
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
             imwrite(imind,cm,'Compare_MA_FDM_homo_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'Compare_MA_FDM_homo_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



% inhomo 1

Length = 10;  % Stablaenge
Time = 5;   % Zetiraum

figure
for n = 1 : 51
    t_fdm = Time / 50 / 0.1 * (n - 1) + 1;
    t_modal = Time / 50 / 0.01 * (n - 1) + 1;
    clf
    plot(x_fdm_inhomo_1, f_fdm_inhomo_1(t_fdm,:),'b-','LineWidth',5)
    hold on
    plot(x_modal_inhomo_1, f_modal_inhomo_1(t_modal,:),'r-','LineWidth',5)
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
    T = ylabel('$f$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(Time/50*(n-1)),'$'];
    T = text(8,1.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'Compare_MA_FDM_inhomo_1.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'Compare_MA_FDM_inhomo_1.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end


% inhomo 2
Length = 10;  % Stablaenge
Time = 20;   % Zetiraum

figure
for n = 1 : 51
    t_fdm = Time / 50 / 0.1 * (n - 1) + 1;
    t_modal = Time / 50 / 0.01  * (n - 1) + 1;
    clf
    plot(x_fdm_inhomo_2, f_fdm_inhomo_2(t_fdm,:),'b-','LineWidth',5)
    hold on
    plot(x_modal_inhomo_2, f_modal_inhomo_2(t_modal,:),'r-','LineWidth',5)
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
    T = ylabel('$f$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(Time/50*(n-1)),'$'];
    T = text(8,1.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'Compare_MA_FDM_inhomo_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'Compare_MA_FDM_inhomo_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end
