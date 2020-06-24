clc
clear
close all

printfigure = 0;

load('Messung')

nt = length(t);
nx = length(x);

sampling_number = 20;
sampling_index = zeros(nt,sampling_number);
sampling_index(:,1) = ceil(rand(nt,1)*nx);
for n = 1 : nt
    k = 2;
    while k <= sampling_number
        temp = ceil(rand()*nx);
        if abs(sampling_index(n,:) - temp) ~= 0
            sampling_index(n,k) = temp;
            k = k + 1;
        end
    end
    sampling_index(n,:) = sort(sampling_index(n,:));
end


Psi = zeros(nx,100);
for i = 1 : nx
    for j = 1 : 100
        Psi(i,j) = sin(2*pi*j*x(i));
    end
end

Phi = zeros(nt, sampling_number, nx);
for n = 1 : nt
    for i = 1 : sampling_number
        Phi(n,i,sampling_index(n,i)) = 1;
    end
end

A = diag(ones(100,1)*0.995);
for i = 1 : 99
    A(i,i+1) = 0.005;
end

z = zeros(100,nt);

for n = 2 : nt
    h = squeeze(Phi(n,:,:)) * Psi;
    y_sampling = y(sampling_index(n,:),n);
    zp = A * z(:, n-1);
    
    e = 0.1;
    cvx_begin;
        variable a(100,1);
        minimize(norm(a,1));
        subject to;
            norm(h * (zp + a) - y_sampling) <= e;
    cvx_end;
    z(:, n) = zp + a;
end


y_re = zeros(nx,nt); 
for n = 1 : nt
    y_re(:,n) =  Psi * z(:,n);
end

figure
for n = 1 : nt
    clf
    plot(x, y(:,n),'k-','LineWidth',3)
    hold on
    plot(x, y_re(:,n),'r-','LineWidth',2)
    hold on
    plot(x(sampling_index(n,:)),y(sampling_index(n,:),n),'c.','Markersize',40)
    ylim([-25 25])
    legend('Signal','recovered Signal','Measure Points')
    set(gca,'Fontsize',20)
    set(gca,'fontname','times new Roman')
    T = title('Signal Recovery','fontsize',40);
    set(T,'Interpreter','latex')
    T = xlabel('$x$','fontsize',30);
    set(T,'Interpreter','latex')
    T = ylabel('$T$','fontsize',30);
    set(T,'Interpreter','latex')
    txt = ['$N = ',num2str(sampling_number),'$'];
    T = text(0.5,22,txt,'FontSize',30);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    name = ['KFCS_TF_2.gif'];
    if printfigure == 1
        if n==1
             imwrite(imind,cm,name,'gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,name,'gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end



