% wavelet approximation

clc
clear
close all

printfigure = 0;

load('Messwerte')

x = 0 : 0.01 : 10-0.01;
nx = length(x);
p = 0 : 0.1 : 10 - 0.1;
m = f(:,1:10:1000);

nt = size(m,1);

number = 6;
S = zeros(nt,number);
S(:,1) = ceil(rand(nt,1)*100);
for t = 1 : nt
    k = 2;
    while k <= number
        temp = ceil(rand()*100);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end


WT = DWT(1000,'haar');
THETA = WT^-1;


Phi = zeros(nt, number, nx);
for t = 1 : nt
    for i = 1 : number
        Phi(t,i,(S(t,i)-1)*10+1) = 1;
    end
end

z = zeros(nx, nt);
z(:,1) = WT * f(1,1:end-1)';

for t = 2 : nt
    h = squeeze(Phi(t,:,:)) * THETA;
    y = m(t,S(t,:))';
    
    e = 0;
    cvx_begin;
        variable a(nx,1);
        minimize(norm(a,2));
        subject to;
            norm(h * (z(:,t-1) + a) - y) <= e;
    cvx_end;
    z(:, t) = z(:,t-1) + a;
end

f_e_wt = zeros(nx,nt);
for n = 1 : nt
    f_e_wt(:,n) = THETA * z(:,n);
end

figure
for n = 1 : 40 : nt
    clf
    plot(x, f(n,1:end-1),'k-','LineWidth',5)
    hold on
    plot(x, f_e_wt(:,n),'c-','LineWidth',5)
    hold on
    plot(x((S(n,:)-1)*10+1),m(n,S(n,:))','r.','Markersize',40)
    legend('TV Real','TV Estimated','Measure Points')
    xlim([0 10])
    ylim([-0.5 2.5])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    z = title('Temperature Distribution','fontsize',40);
    set(z,'Interpreter','latex')
    z = xlabel('$x$','fontsize',30);
    set(z,'Interpreter','latex')
    z = ylabel('$T$','fontsize',30);
    set(z,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    z = text(0.8,0.6,txt,'FontSize',30);
    set(z,'Interpreter','latex')
    txt = ['$N = ',num2str(number),'$'];
    z = text(0.8,0.8,txt,'FontSize',30);
    set(z,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['KFCS_wavelet_6_l2.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

save('kfcs_tv_wavelet_6_l2.mat','f_e_wt','S')




