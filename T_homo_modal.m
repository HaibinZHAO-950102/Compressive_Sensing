clc
clear
printfigure = 0;

Length = 1;  % length of the rod
Time = 0.05;  % duration
step_length = 0.001;
step_time = 0.0001;
x = 0 : step_length : Length;
t = 0 : step_time : Time;
N_length = length(x);
N_time = length(t);

k = 10;  % thermal conductivity in cm^2/100s

N = 10;  % order of eigen functions
lambda = 0 : pi / Length : N * pi / Length;  % frequence of eigen functions
f = zeros(N_time, N_length);  % temperature
f(1,:) = sin(x / Length * 2 * pi);
% f(1,:) = floor(x * 5) / 5;
phi = zeros(N + 1, N_length);  % eigen function
T = zeros(N + 1, N_time);  % weights

step_time_max = 2 / k / lambda(end)^2;
if step_time > step_time_max
    ['step_time should be under ',num2str(step_time_max)]
end


phi(1,:) = sqrt(1 / Length);
for i = 2 : N + 1
    for n = 1 : N_length
        phi(i, n) = sqrt(2 / Length) * cos(lambda(i) * x(n));
    end
end

for i = 1 : N + 1
    plot(0:step_length:Length,phi(i,:))
    hold on
end
setplt('Eigenfunctions','$x$','$value$','Eigenfunktionen',printfigure)

for i = 1 : N + 1
    T(i, 1) = 0;
    for n = 1 : N_length
        T(i, 1) = T(i, 1) + f(1, n) * phi(i, n) * step_length;
    end
end
for i = 1 : N + 1
    for n = 2 : N_time
        T(i, n) = (1 - step_time * k * lambda(i)^2)^n * T(i, 1);
    end
end
for n = 2 : N_time
    for i = 1 : N + 1
        f(n,:) = f(n,:) + T(i,n) * phi(i,:);
    end
end

figure
for n = 1 : 5 : N_time
    plot(x, f(n,:))
    ylim([-1 1])
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Temperature Distribution','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    T = text(0.8,0.6,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = ',num2str(N),'$'];
    T = text(0.8,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if n==1
         imwrite(imind,cm,'Temperature_Distribution.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
    else
         imwrite(imind,cm,'Temperature_Distribution.gif','gif','WriteMode','append','DelayTime',1e-6);
    end
end

figure
[X, Y] = meshgrid(x, t);
mesh(X,Y,f)
setmesh('Tempreature Distribution','$x$','$t$','$T$','T-t-x',printfigure)








