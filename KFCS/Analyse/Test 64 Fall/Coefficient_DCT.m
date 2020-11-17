clc
clear
close all
printfigure = 0;

load('Messwerte.mat')

Dt = 0.1;

y = f(p_index,1:Dt/dt:end);
x = 0 : dx : 10;
x = x(p_index);
nt = size(y,2);
nx = length(x);

% WT = DWT(nx, 'db4');
DCT = idct(eye(nx))^-1;
Psi = DCT^-1;

z = zeros(nt,nx);

for i = 1 : nt
    z(i,:) = (DCT * y(:,i))';
end

figure
for n = 1 : 0.5/Dt : nt
    clf
    plot(z(n,:),'k.','Markersize',20)
    xlim([1 nx])
    ylim([min(min(z)) max(max(z))])
    setplt('Coefficient WT','$z$','$value$','Temperature Distribution',0)
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(800, 20,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['Koeffizienten_WT.gif'];
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
for n = 1 : 0.5/Dt : nt-1
    clf
    plot(dz(n,:),'k.','Markersize',20)
    xlim([1 nx])
    ylim([min(min(dz)) max(max(dz))])
    setplt('Sorted Coefficient Change','$\Delta z$','$value$','Temperature Distribution',0)
    txt = ['$t = ',num2str((n-1)*Dt),'$'];
    T = text(800,0.04,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['Koeffizienten_WT_Change.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



