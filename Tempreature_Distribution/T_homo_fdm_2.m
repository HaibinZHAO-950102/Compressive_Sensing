clc
clear
close all

printfigure = 0;

Length = 10;  % Stablaenge
Time = 0.5;   % Zetiraum

step_length = 0.1;
step_time = 0.01;

x = 0 : step_length : Length;
t = 0 : step_time : Time;

N_length = length(x);
N_time = length(t);

k = 0.1;

f = zeros(N_time, N_length);
f(1,:) = floor(x / 2) / 5;

A = step_length ^ 2 / (step_length ^ 2 + 2 * k * step_time);
B = (k * step_time) / (step_length ^ 2 + 2 * k * step_time);

H = speye(N_length * N_time);

for n = 2 : N_length - 1
    for m = 2 : N_time
        H(n+N_length*(m-1),n+N_length*(m-1)) = -1;
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

b = zeros(N_length*N_time,1);
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
for n = 1 : 1 : N_time
    plot(x, f(n,:),'LineWidth',5)
    ylim([0 1])
    xlim([0 10])
    setplt('Temperature Distribution','$x$','$f$','Temperature Distribution',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    T = text(1,0.8,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'TV_homo_fdm_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'TV_homo_fdm_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

figure
[X, Y] = meshgrid(x, t);
mesh(X,Y,f)
setmesh('Tempreature Distribution','$x$','$t$','$f$','TV_homo_fdm_Ttx_2',printfigure)
    

f_fdm_homo_2 = f;
x_fdm_homo_2 = x;
save('f_fdm_homo_2.mat','x_fdm_homo_2','f_fdm_homo_2')
