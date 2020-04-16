clc
clear
printfigure = 1;

Length = 1;  % length of the rod
Time = 5;  % duration
step_length = 0.001;
step_time = 0.01;
x = 0 : step_length : Length;
t = 0 : step_time : Time;
N_length = length(x);
N_time = length(t);

k = 0.1;  % thermal conductivity in cm^2/s


N = 10;  % order of eigen functions
lambda = 0 : pi / Length : N * pi / Length;  % frequence of eigen functions
f = zeros(N_time, N_length);  % temperature
f(1,:) = sin(x / Length * 2 * pi) + 1;
% f(1,:) = floor(x * 5) / 5;
phi = zeros(N + 1, N_length);  % eigen function
T = zeros(N + 1, N_time);  % weights
u = zeros(N_time, N_length);
U = zeros(N + 1, N_time);  % weights

plot(x,f(1,:))
setplt('Initial Condition','$x$','$T$','TV_inhomo_modal_inital_condition',0)

phi(1,:) = sqrt(1 / Length);
for i = 2 : N + 1
    for n = 1 : N_length
        phi(i, n) = sqrt(2 / Length) * cos(lambda(i) * x(n));
    end
end

figure
for i = 1 : N + 1
    plot(0:step_length:Length,phi(i,:))
    hold on
end
txt = ['$N = ',num2str(N),'$'];
TEXT = text(0.8,1.2,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setplt('Eigenfunctions','$x$','$value$','TV_inhomo_modal_Eigenfunctions',0)

[index,minimum] = min(abs(x - 0.3));
u(:,301) = 50 * t;
for i = 1 : N + 1
    for j = 1 : N_time
        U(i, j) = 0;
        for n = 1 : N_length
            U(i, j) = U(i, j) + u(j, n) * phi(i, n) * step_length;
        end
    end
end

for i = 1 : N + 1
    T(i, 1) = 0;
    for n = 1 : N_length
        T(i, 1) = T(i, 1) + f(1, n) * phi(i, n) * step_length;
    end
end

for i = 1 : N + 1
    for n = 2 : N_time
        T(i, n) = (1 - step_time * k * lambda(i)^2) * T(i, n - 1) + step_time * U(i, n);
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
    ylim([0 2])
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
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'TV_inhomo_modal_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'TV_inhomo_modal_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

figure
[X, Y] = meshgrid(x, t);
mesh(X,Y,f)
setmesh('Tempreature Distribution','$x$','$t$','$T$','TV_inhomo_modal_Ttx_2',printfigure)








