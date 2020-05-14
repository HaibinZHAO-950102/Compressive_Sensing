clc
clear
close all
printfigure = 1;

A = diag(ones(100,1)*0.995);
for i = 1 : 99
    A(i,i+1) = 0.005;
end

t = 0 : 0.01 : 0.5;
x = 0 : 0.0001 : 1;

nt = length(t);
nx = length(x);

z = zeros(100,nt);
z(85,1) = 10;
z(25,1) = 10;

u = zeros(100,nt);
u(100,:) = 2*sin(2*pi*nt/0.5);
u(50,:) = -0.1*cos(2*pi*2*nt/0.5);


for i = 2 : nt
    z(:,i) = A * z(:,i-1) + u(:,i-1);
end

Psi = zeros(nx,100);
for i = 1 : nx
    for j = 1 : 100
        Psi(i,j) = sin(2*pi*j*x(i));
    end
end

y = zeros(nx,nt);
for i = 1 : nt
    y(:,i) = Psi * z(:,i);
end

for n = 1 : nt
    clf
    plot(x, y(:,n),'k-','LineWidth',1)
    ylim([-25 25])
    setplt('Test Function','$x$','$y$','Test Function',0)
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'Test Function.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'Test Function.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

save('Messung','y','x','t')

        