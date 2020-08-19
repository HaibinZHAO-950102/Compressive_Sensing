clc
clear
close all
printfigure = 0;

load Messwerte
load f_e_scr_kf_fdm
load f_e_scr_kf_fdm_ndw

Dt = 0.1;

f = f(:,1:Dt/dt:end);
f1 = f_e_scr_kf_fdm_ndw;
f2 = f_e_scr_kf_fdm;


NAME = 'comapre dynamic weight';
l1 = 'KFCS';
l2 = 'KFCS dynamic weight';




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

for i = 1 : size(f,2)
    similariy_1(i) = f(:,i)' * f1(:,i) / (norm(f(:,i)) * norm(f1(:,i)));
    similariy_2(i) = f(:,i)' * f2(:,i) / (norm(f(:,i)) * norm(f2(:,i)));
end

error_1 = 2 * (1-similariy_1);
error_2 = 2 * (1-similariy_2);

mean_error_1 = mean(error_1(25:end));
mean_error_2 = mean(error_2(25:end));

t = 0 : Dt : 20;
figure
plot(t(25:end),error_1(25:end),'r-','LineWidth',5)
hold on
plot(t(25:end),error_2(25:end),'b-','LineWidth',5)
hold on
plot(t(25:end),ones(length(t(25:end)))*mean_error_1,'r-.','LineWidth',2)
hold on
plot(t(25:end),ones(length(t(25:end)))*mean_error_2,'b-.','LineWidth',2)
hold on
legend(l1,l2)
setplt('Error','$t$','Error',[NAME,' error'],printfigure)
