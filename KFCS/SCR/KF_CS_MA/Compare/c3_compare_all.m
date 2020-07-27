clc
clear
close all
printfigure = 0;

load f_random_simpling_KF_modal
load f_e_mt_kf_modal
load f_e_mt_kf_modal_integrated_1
load f_e_kalman_modal_1D_12
load f_e_kalman_modal_1D_20
load Messwerte

Dt = 0.1;

f = f(:,1:Dt/dt:end);
f2 = f_random_simpling_KF_modal;
f5 = f_e_mt_kf_modal_integrated_1;
f1 = f_e_kalman_modal_1D_12;
f4 = f_e_kalman_modal_1D_20;
f3 = f_e_mt_kf_modal;


NAME = 'comapre all';
l2 = 'random sampling KF';
l5 = 'KFCS integrated';
l1 = 'KF 12';
l4 = 'KF 20';
l3 = 'KFCS';





x = 0 : dx : 10;

% figure
% for n = 1 : 0.5/Dt : size(f,2)
%     clf
%     plot(x, f(:,n),'k-','LineWidth',5)
%     hold on
%     plot(x, f1(:,n),'r-','LineWidth',5)
%     hold on
%     plot(x, f2(:,n),'b-','LineWidth',5)
%     legend('Origin',l1,l2)
%     xlim([0 10])
%     ylim([-0.5 2.5])
%     set(gca,'Fontsize',20)
%     set(gca,'fontname','times new Roman')
%     T = title('Temperature Distribution','fontsize',40);
%     set(T,'Interpreter','latex')
%     T = xlabel('$x$','fontsize',30);
%     set(T,'Interpreter','latex')
%     T = ylabel('$f$','fontsize',30);
%     set(T,'Interpreter','latex')
%     set(gcf,'outerposition',get(0,'screensize'));
%     txt = ['$t = ',num2str(Dt * (n-1)),'$'];
%     T = text(0.8,0.6,txt,'FontSize',30);
%     set(T,'Interpreter','latex')
%     txt = ['$M = 6$'];
%     T = text(0.8,0.8,txt,'FontSize',30);
%     set(T,'Interpreter','latex')
%     drawnow
%     frame=getframe(gcf);
%     imind=frame2im(frame);
%     [imind,cm] = rgb2ind(imind,256);
%     name = [NAME,'.gif'];
%     if printfigure == 1
%         if n==1
%              imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
%         else
%              imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
%         end
%     end
% end



error_1 = f - f1;
error_2 = f - f2;
error_3 = f - f3;
error_4 = f - f4;
error_5 = f - f5;

zero = zeros(size(f));
% figure
% for n = 1 : 0.5/Dt : size(f,2)
%     clf
%     plot(x, error_1(:,n),'r-','LineWidth',5)
%     hold on
%     plot(x, error_2(:,n),'b-','LineWidth',5)
%     plot(x, zero(:,n),'k--','LineWidth',3)
%     legend(l1,l2)
%     xlim([0 10])
%     ylim([-1 1])
%     set(gca,'Fontsize',20)
%     set(gca,'fontname','times new Roman')
%     T = title('Error','fontsize',40);
%     set(T,'Interpreter','latex')
%     T = xlabel('$x$','fontsize',30);
%     set(T,'Interpreter','latex')
%     T = ylabel('$f$','fontsize',30);
%     set(T,'Interpreter','latex')
%     set(gcf,'outerposition',get(0,'screensize'));
%     txt = ['$t = ',num2str(Dt * (n-1)),'$'];
%     T = text(0.8,0.6,txt,'FontSize',30);
%     set(T,'Interpreter','latex')
%     txt = ['$M = 12$'];
%     T = text(0.8,0.8,txt,'FontSize',30);
%     set(T,'Interpreter','latex')
%     drawnow
%     frame=getframe(gcf);
%     imind=frame2im(frame);
%     [imind,cm] = rgb2ind(imind,256);
%     name = [NAME,' error.gif'];
%     if printfigure == 1
%         if n==1
%              imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
%         else
%              imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
%         end
%     end
% end

similariy_1 = zeros(size(f,2),1);
similariy_2 = zeros(size(f,2),1);
similariy_3 = zeros(size(f,2),1);
similariy_4 = zeros(size(f,2),1);
similariy_5 = zeros(size(f,2),1);


for i = 1 : size(f,2)
    similariy_1(i) = f(:,i)' * f1(:,i) / (norm(f(:,i)) * norm(f1(:,i)));
    similariy_2(i) = f(:,i)' * f2(:,i) / (norm(f(:,i)) * norm(f2(:,i)));
    similariy_3(i) = f(:,i)' * f3(:,i) / (norm(f(:,i)) * norm(f3(:,i)));
    similariy_4(i) = f(:,i)' * f4(:,i) / (norm(f(:,i)) * norm(f4(:,i)));
    similariy_5(i) = f(:,i)' * f5(:,i) / (norm(f(:,i)) * norm(f5(:,i)));
end

error_1 = 2 * (1-similariy_1);
error_2 = 2 * (1-similariy_2);
error_3 = 2 * (1-similariy_3);
error_4 = 2 * (1-similariy_4);
error_5 = 2 * (1-similariy_5);

mean_error_1 = mean(error_1(10:end));
mean_error_2 = mean(error_2(10:end));
mean_error_3 = mean(error_3(10:end));
mean_error_4 = mean(error_4(10:end));
mean_error_5 = mean(error_5(10:end));

t = 0 : Dt : 20;
figure
plot(t(10:end),error_1(10:end),'r-','LineWidth',5)
hold on
plot(t(10:end),error_2(10:end),'b-','LineWidth',5)
hold on
plot(t(10:end),error_3(10:end),'g-','LineWidth',5)
hold on
plot(t(10:end),error_4(10:end),'c-','LineWidth',5)
hold on
plot(t(10:end),error_5(10:end),'m-','LineWidth',5)
hold on
plot(t(10:end),ones(length(t(10:end)))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t(10:end),ones(length(t(10:end)))*mean_error_2,'b-.','LineWidth',2)
hold on
plot(t(10:end),ones(length(t(10:end)))*mean_error_3,'g-.','LineWidth',2)
hold on
plot(t(10:end),ones(length(t(10:end)))*mean_error_4,'c-.','LineWidth',2)
hold on
plot(t(10:end),ones(length(t(10:end)))*mean_error_5,'m-.','LineWidth',2)
hold on
legend(l1,l2,l3,l4,l5)
setplt('Error','$t$','Error',[NAME,' error'],printfigure)
