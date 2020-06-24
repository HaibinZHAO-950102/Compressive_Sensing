clc
clear
close all
printfigure = 1;

load('Messwerte')
load('f_e_wt_kf_naiv_1')
load('f_e_wt_kf_naiv_50')
load('f_e_wt_kf_naiv_100')
load('f_e_wt_kf_naiv_500')
load('f_e_wt_kf_naiv_1000')
load('f_e_wt_kf_naiv_5000')

f = f(:,1:end-1);

x = 0 : 0.01 :10 - 0.01;

figure
for n = 1 : (size(f,1)-1)/50 : (size(f,1)-1)
    clf
    plot(x, f(n,:),'k-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf_naiv_1(:,n),'r-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf_naiv_50(:,n),'g-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf_naiv_100(:,n),'b-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf_naiv_500(:,n),'c-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf_naiv_1000(:,n),'m-','LineWidth',3)
    hold on
    plot(x, f_e_wt_kf_naiv_5000(:,n),'y-','LineWidth',3)
    hold on
    legend('Origin','G1/G2 1','G1/G2 50','G1/G2 100','G1/G2 500','G1/G2 1000','G1/G2 5000')
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
    txt = ['$t = ',num2str((n-1)/(size(f,1)-1)*20),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = 6$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['compare KFCS_Naiv.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



similariy_1 = zeros(size(f,1),1);
similariy_50 = zeros(size(f,1),1);
similariy_100 = zeros(size(f,1),1);
similariy_500 = zeros(size(f,1),1);
similariy_1000 = zeros(size(f,1),1);
similariy_5000 = zeros(size(f,1),1);

for i = 1 : size(f,1)
    similariy_1(i) = f(i,:) * f_e_wt_kf_naiv_1(:,i) / (norm(f(i,:)) * norm(f_e_wt_kf_naiv_1(:,i)));
    similariy_50(i) = f(i,:) * f_e_wt_kf_naiv_50(:,i) / (norm(f(i,:)) * norm(f_e_wt_kf_naiv_50(:,i)));
    similariy_100(i) = f(i,:) * f_e_wt_kf_naiv_100(:,i) / (norm(f(i,:)) * norm(f_e_wt_kf_naiv_100(:,i)));
    similariy_500(i) = f(i,:) * f_e_wt_kf_naiv_500(:,i) / (norm(f(i,:)) * norm(f_e_wt_kf_naiv_500(:,i)));
    similariy_1000(i) = f(i,:) * f_e_wt_kf_naiv_1000(:,i) / (norm(f(i,:)) * norm(f_e_wt_kf_naiv_1000(:,i)));
    similariy_5000(i) = f(i,:) * f_e_wt_kf_naiv_5000(:,i) / (norm(f(i,:)) * norm(f_e_wt_kf_naiv_5000(:,i)));
end

error_1 = 2 * (1-similariy_1);
error_50 = 2 * (1-similariy_50);
error_100 = 2 * (1-similariy_100);
error_500 = 2 * (1-similariy_500);
error_1000 = 2 * (1-similariy_1000);
error_5000 = 2 * (1-similariy_5000);

t = 0 : 20/(size(f,1)-1) : 20;
figure
plot(t(10:end),error_1(10:end),'r-','LineWidth',3)
hold on
plot(t(10:end),error_50(10:end),'g-','LineWidth',3)
hold on
plot(t(10:end),error_100(10:end),'b-','LineWidth',3)
hold on
plot(t(10:end),error_500(10:end),'c-','LineWidth',3)
hold on
plot(t(10:end),error_1000(10:end),'m-','LineWidth',3)
hold on
plot(t(10:end),error_5000(10:end),'y-','LineWidth',3)
hold on
ylim([0 0.02])
legend('G1/G2 1','G1/G2 50','G1/G2 100','G1/G2 500','G1/G2 1000','G1/G2 5000')
setplt('Error','$t$','Error','KFCS_Naiv_Error',printfigure)
