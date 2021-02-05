clc
clear
close all
printfigure = 1;

load Messwerte
load f_e_wt_kf_naiv_kkk1
load f_e_wt_kf_naiv_kkk10
load f_e_wt_kf_naiv_kkk100
load f_e_wt_kf_naiv_kkk1000
load f_e_wt_kf_naiv_kkk10000
load f_e_wt_kf_naiv_kkk100000
load f_e_wt_kf_naiv_kkk1000000
load f_e_wt_kf_naiv_kkk10000000
load f_e_wt_kf_naiv_kkk100000000
load f_e_wt_kf_naiv_kkk1000000000

Dt = 0.1;

f = f(:,1:Dt/dt:end);
f1 = f_e_wt_kf_naiv_kkk1;
f2 = f_e_wt_kf_naiv_kkk10;
f3 = f_e_wt_kf_naiv_kkk100;
f4 = f_e_wt_kf_naiv_kkk1000;
f5 = f_e_wt_kf_naiv_kkk10000;
f6 = f_e_wt_kf_naiv_kkk100000;
f7 = f_e_wt_kf_naiv_kkk1000000;
f8 = f_e_wt_kf_naiv_kkk10000000;
f9 = f_e_wt_kf_naiv_kkk100000000;
f10 = f_e_wt_kf_naiv_kkk1000000000;

NAME = 'comapre g2g1';
l1 = 'G2/G1 = 10^0';
l2 = 'G2/G1 = 10^1';
l3 = 'G2/G1 = 10^2';
l4 = 'G2/G1 = 10^3';
l5 = 'G2/G1 = 10^4';
l6 = 'G2/G1 = 10^5';
l7 = 'G2/G1 = 10^6';
l8 = 'G2/G1 = 10^7';
l9 = 'G2/G1 = 10^8';
l10 = 'G2/G1 = 10^9';




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
error_6 = f - f6;
error_7 = f - f7;
error_8 = f - f8;
error_9 = f - f9;
error_10 = f - f10;


% zero = zeros(size(f));
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
similariy_6 = zeros(size(f,2),1);
similariy_7 = zeros(size(f,2),1);
similariy_8 = zeros(size(f,2),1);
similariy_9 = zeros(size(f,2),1);
similariy_10 = zeros(size(f,2),1);


for i = 1 : size(f,2)
    similariy_1(i) = f(:,i)' * f1(:,i) / (norm(f(:,i)) * norm(f1(:,i)));
    similariy_2(i) = f(:,i)' * f2(:,i) / (norm(f(:,i)) * norm(f2(:,i)));
    similariy_3(i) = f(:,i)' * f3(:,i) / (norm(f(:,i)) * norm(f3(:,i)));
    similariy_4(i) = f(:,i)' * f4(:,i) / (norm(f(:,i)) * norm(f4(:,i)));
    similariy_5(i) = f(:,i)' * f5(:,i) / (norm(f(:,i)) * norm(f5(:,i)));
    similariy_6(i) = f(:,i)' * f6(:,i) / (norm(f(:,i)) * norm(f6(:,i)));
    similariy_7(i) = f(:,i)' * f7(:,i) / (norm(f(:,i)) * norm(f7(:,i)));
    similariy_8(i) = f(:,i)' * f8(:,i) / (norm(f(:,i)) * norm(f8(:,i)));
    similariy_9(i) = f(:,i)' * f9(:,i) / (norm(f(:,i)) * norm(f9(:,i)));
    similariy_10(i) = f(:,i)' * f10(:,i) / (norm(f(:,i)) * norm(f10(:,i)));
end

error_1 = 2 * (1-similariy_1);
error_2 = 2 * (1-similariy_2);
error_3 = 2 * (1-similariy_3);
error_4 = 2 * (1-similariy_4);
error_5 = 2 * (1-similariy_5);
error_6 = 2 * (1-similariy_6);
error_7 = 2 * (1-similariy_7);
error_8 = 2 * (1-similariy_8);
error_9 = 2 * (1-similariy_9);
error_10 = 2 * (1-similariy_10);

mean_error_1 = mean(error_1);
mean_error_2 = mean(error_2);
mean_error_3 = mean(error_3);
mean_error_4 = mean(error_4);
mean_error_5 = mean(error_5);
mean_error_6 = mean(error_6);
mean_error_7 = mean(error_7);
mean_error_8 = mean(error_8);
mean_error_9 = mean(error_9);
mean_error_10 = mean(error_10);

t = 0 : Dt : 20;
figure
figure1 = plot(t,error_1,'r-','LineWidth',5)
hold on
figure3 = plot(t,error_3,'b-','LineWidth',5)
hold on
figure5 = plot(t,error_5,'b-','LineWidth',5)
hold on
figure7 = plot(t,error_7,'b-','LineWidth',5)
hold on
figure9 = plot(t,error_9,'b-','LineWidth',5)
hold on
figure11 = plot(t,ones(length(t))*mean_error_1,'r-.','LineWidth',2)
hold on
figure13 = plot(t,ones(length(t))*mean_error_3,'b-.','LineWidth',2)
hold on
figure15 = plot(t,ones(length(t))*mean_error_5,'b-.','LineWidth',2)
hold on
figure17 = plot(t,ones(length(t))*mean_error_7,'b-.','LineWidth',2)
hold on
figure19 = plot(t,ones(length(t))*mean_error_9,'b-.','LineWidth',2)
hold on
set(figure1,'color',[255,0,0]/255);
set(figure3,'color',[255,255,0]/255);
set(figure5,'color',[0,0,0]/255);
set(figure7,'color',[0,0,255]/255);
set(figure9,'color',[255,0,255]/255);
set(figure11,'color',[255,0,0]/255);
set(figure13,'color',[255,255,0]/255);
set(figure15,'color',[0,0,0]/255);
set(figure17,'color',[0,0,255]/255);
set(figure19,'color',[255,0,255]/255);

legend(l1,l3,l5,l7,l9)
ylim([0 0.005])
setplt('','$t$','$\epsilon$',[NAME,' error'],printfigure)

MEANERROR = [mean_error_1,mean_error_2,mean_error_3,mean_error_4...
    ,mean_error_5,mean_error_6,mean_error_7,mean_error_8,mean_error_9...
    ,mean_error_10];

figure
plot(0:9,MEANERROR,'b-','linewidth',3)
hold on
plot(0:9,MEANERROR,'r.','markersize',40)
xlim([0 9])
setplt('','$log_{10}{G2/G1}$','$\bar{\epsilon}$','mean_error_g2g1',printfigure)
    xticks(0:1:9)
    xticklabels({'0','1','2','3','4','5','6','7','8','9'})

close all