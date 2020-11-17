clc
clear
close all

printfigure = 1;

load('Messwerte')


Dt = 0.1;
m = m(:,1:Dt/dt:end);

N_time = size(m,2);

number = 12;
S = zeros(N_time,number);
S(:,1) = ceil(rand(N_time,1)*64);
for t = 1 : N_time
    k = 2;
    while k <= number
        temp = ceil(rand()*64);
        if abs(S(t,:) - temp) ~= 0
            S(t,k) = temp;
            k = k + 1;
        end
    end
    S(t,:) = sort(S(t,:));
end

M = number;

x = 0 : dx : 10;
nx = length(x);
t = 0 : Dt : 20;
nt = length(t);


t = 0 : 0.1:20;
for n = 0 : 10
    figure
    timef = 2000 / 10 * n + 1
    timefe = 200 / 10 * n + 1
    plot(x(p_index(S(timefe,:))),m(S(timefe,:),timefe),'r.','Markersize',40)
    legend('Messungen')

    xlim([0 10])
    ylim([-0.5 2.5])
    txt = ['$t = ',num2str(2*n),'$'];
    T = text(0.8,0.4,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['RS_M_shot_',num2str(n+1)];
    setplt('Temperaturverteilung','$x$','$f$',name,printfigure)
end

for n = 0 : 10
    figure
    timef = 2000 / 10 * n + 1
    timefe = 200 / 10 * n + 1
    tv = plot(x, f(:,timef),'k-','LineWidth',5);
    set(tv,'color',[0.8,0.8,0.8])
    hold on
    plot(x(p_index(S(timefe,:))),m(S(timefe,:),timefe),'r.','Markersize',40)
    legend('Signal','Messungen')

    xlim([0 10])
    ylim([-0.5 2.5])
    txt = ['$t = ',num2str(2*n),'$'];
    T = text(0.8,0.4,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    txt = ['$M = ',num2str(number),'$'];
    T = text(0.8,0.8,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['RS_shot_',num2str(n+1)];
    setplt('Temperaturverteilung','$x$','$f$',name,printfigure)
end

close all
