clc
clear
close all

printfigure = 1;

a = 0.1;

Length = 10;  % Stablaenge
Time = 5;   % Zetiraum

step_length = 0.1;
step_time = 0.1;

x = 0 : step_length : Length;
t = 0 : step_time : Time;

N_length = length(x);
N_time = length(t);

k = 0.1;

f = zeros(N_time, N_length);
f(1,:) = sin(x / Length * 2 * pi) + 1;

u = zeros(N_time, N_length);

A = step_length ^ 2 / (step_length ^ 2 + 2 * k * step_time);
B = (k * step_time) / (step_length ^ 2 + 2 * k * step_time);
C =  - step_length ^ 2  * step_time / (step_length ^ 2 + 2 * k * step_time);
K = a * C - 1;

H = speye(N_length * N_time);
b = zeros(N_length*N_time,1);

for n = 2 : N_length - 1
    for m = 2 : N_time
        H(n+N_length*(m-1),n+N_length*(m-1)) = K;
        H(n+N_length*(m-1),n+N_length*(m-2)) = A;
        H(n+N_length*(m-1),n+N_length*(m-1)-1) = B;
        H(n+N_length*(m-1),n+N_length*(m-1)+1) = B;
    end
end

length_edge = [1, N_length];
for n = 1 : 2
    for m = 2 : N_time
        H(length_edge(n)+N_length*(m-1),length_edge(n)+N_length*(m-1)) = -1;
        H(length_edge(n)+N_length*(m-1),length_edge(n)+N_length*(m-1)-sign(length_edge(n)-N_length/2)) = 1;
    end
end

for n = 1 : N_length
    H(n,n) = 1;
    b(n) = f(1,n);
end

T = H\b;
for n = 1 : N_length
    for m = 2 : N_time
        f(m,n) = T(n+N_length*(m-1));
    end
end

figure
for n = 1 : N_time
    plot(x, f(n,:),'LineWidth',5)
    ylim([0 2])
    setplt('','$x$','$f$','Temperature Distribution',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    T = text(8,1.6,txt,'FontSize',60);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'TV_inhomo_fdm_1.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'TV_inhomo_fdm_1.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

figure
[X, Y] = meshgrid(x, t);
mesh(X,Y,f)
setmesh('','$x$','$t$','$f$','TV_inhomo_fdm_1_Ttx',printfigure)

% f_fdm_inhomo_1 = f;
% x_fdm_inhomo_1 = x;
% save('f_fdm_inhomo_1.mat','x_fdm_inhomo_1','f_fdm_inhomo_1')

close all




