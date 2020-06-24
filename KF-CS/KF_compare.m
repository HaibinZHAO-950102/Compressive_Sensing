clc
clear
close all
printfigure = 0;

load('KFCS')
load('f_e_kalman_6')
load('f_e_kalman_11')
load('Messwerte')
load('f_e_wt_kf')

f = f(:,1:end-1);
f_e = f_e(1:end-1,:);
f_e_kalman_6 = f_e_kalman_6(:,1:end-1);
f_e_kalman_11 = f_e_kalman_11(:,1:end-1);

x = 0 : 0.01 :10 - 0.01;

figure
for n = 1 : (size(f_e,2)-1)/50 : size(f_e,2)
    clf
    plot(x, f(n,:),'k-','LineWidth',3)
    hold on
    plot(x, f_e(:,n),'r-','LineWidth',3)
    hold on
    plot(x, f_e_kalman_6(n,:),'b-','LineWidth',3)
    hold on
    plot(x, f_e_kalman_11(n,:),'c-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf(:,n),'Color',[0.1961 0.6314 0.5373],'LineWidth',3)
    hold on
    legend('Origin','KF-CS 6 Sensors','KF 6 Sensors','KF 11 Sensors','CS WT KF 6')
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
    txt = ['$t = ',num2str((n-1)/(size(f_e,2)-1)*20),'$'];
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



error_kfcs = f - f_e';
error_kf_6 = f - f_e_kalman_6;
error_kf_11 = f - f_e_kalman_11;
error_wt_kf = f - f_e_wt_kf';

zero = zeros(size(f));
figure
for n = 1 : (size(f_e,2)-1)/50 : size(f_e,2)
    clf
    plot(x, error_kfcs(n,:),'r-','LineWidth',3)
    hold on
    plot(x, error_kf_6(n,:),'b-','LineWidth',3)
    hold on
    plot(x, error_kf_11(n,:),'c-','LineWidth',3)
    hold on
    plot(x, error_wt_kf(n,:),'Color',[0.1961 0.6314 0.5373],'LineWidth',3)
    hold on
    plot(x, zero(n,:),'k--','LineWidth',5)
    legend('KF-CS','KF 6','KF 11','CS WT KF')
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
    txt = ['$t = ',num2str((n-1)/(size(f_e,2)-1)*20),'$'];
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

similariy_kfcs = zeros(size(f_e,2),1);
similariy_kf_6 = zeros(size(f_e,2),1);
similariy_kf_11 = zeros(size(f_e,2),1);
similariy_wt_kf = zeros(size(f_e,2),1);

for i = 1 : size(f_e,2)
    similariy_kfcs(i) = f(i,:) * f_e(:,i) / (norm(f(i,:)) * norm(f_e(:,i)));
    similariy_kf_6(i) = f(i,:) * f_e_kalman_6(i,:)' / (norm(f(i,:)) * norm(f_e_kalman_6(i,:)));
    similariy_kf_11(i) = f(i,:) * f_e_kalman_11(i,:)' / (norm(f(i,:)) * norm(f_e_kalman_11(i,:)));
    similariy_wt_kf(i) = f(i,:) * f_e_wt_kf(:,i) / (norm(f(i,:)) * norm(f_e_wt_kf(:,i)));
end

error_kfcs = 2 * (1-similariy_kfcs);
error_kf_6 = 2 * (1-similariy_kf_6);
error_kf_11 = 2 * (1-similariy_kf_11);
error_wt_kf = 2 * (1-similariy_wt_kf);

t = 0 : 20/(size(f_e,2)-1) : 20;
figure
plot(t(10:end),error_kfcs(10:end),'r-','LineWidth',3)
hold on
plot(t(10:end),error_kf_6(10:end),'b-','LineWidth',3)
hold on
plot(t(10:end),error_kf_11(10:end),'c-','LineWidth',3)
hold on
plot(t(10:end),error_wt_kf(10:end),'Color',[0.1961 0.6314 0.5373],'LineWidth',3)
legend('KF-CS','KF 6','KF 11','CS WT KF 6')
setplt('Error','$t$','Error','Error',printfigure)
