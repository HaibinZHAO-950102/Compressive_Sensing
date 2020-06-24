clc
clear
close all
printfigure = 1;

load('Messwerte.mat')

y = f(:,1:end-1);
x = 0 : 0.01 : 10-0.001;
nt = size(y,1);
nx = length(x);

WT = DWT(nx, 'haar');
Psi = WT^-1;

z = zeros(nt,nx);

for i = 1 : nt
    z(i,:) = (WT * y(i,:)')';
end

figure
for n = 1 : nt
    clf
    plot(z(n,:),'k.','Markersize',20)
    ylim([-10 30])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Coefficient TV WT','fontsize',40);
    set(T,'Interpreter','latex')
    T = ylabel('$z$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['TV_koeffizienten_WT.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

dz = zeros(nt-1,nx);
for i = 1 : nt-1
    dz(i,:) = sort(abs(z(i+1,:) - z(i,:)),'descend');
end

figure
for n = 1 : nt-1
    clf
    plot(dz(n,:),'k.','Markersize',20)
    ylim([-0.1 0.1])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Sorted Coefficient_change TV WT','fontsize',40);
    set(T,'Interpreter','latex')
    T = ylabel('$dz$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['TV_koeffizienten_WT_change.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



