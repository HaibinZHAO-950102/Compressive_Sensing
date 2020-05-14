clc
clear
close all
printfigure = 0;

load('Messung')
nt = length(t);
nx = length(x);

sampling_number = 401;
M = sampling_number;
% sampling_index = zeros(nt,sampling_number);
% sampling_index(:,1) = ceil(rand(nt,1)*nx);
% 
% for i = 1 : nt
%     k = 2;
%     while k <= sampling_number
%         temp = ceil(rand()*nx);
%         if abs(sampling_index(i,:)-temp) ~= 0
%             sampling_index(i,k) = temp;
%             k = k + 1;
%         end
%     end
%     sampling_index(i,:) = sort(sampling_index(i,:));
% end

sampling_index = zeros(nt,sampling_number);
for i = 1 : nt
    sampling_index(i,:) = 1 : 25 : 10001;
end

y_sampling = zeros(sampling_number,nt);
for i = 1 : nt
    y_sampling(:,i) = y(sampling_index(i,:),i);
end

Psi = zeros(nx,100);
for i = 1 : nx
    for j = 1 : 100
        Psi(i,j) = sin(2*pi*j*x(i));
    end
end

Phi = zeros(sampling_number, nx, nt);
for i = 1 : nt
    for j = 1 : sampling_number
        Phi(j,sampling_index(i,j),i) = 1;
    end
end

N = 100;
z = zeros(N,nt);
z(85,1) = 10;
z(25,1) = 10;

A = diag(ones(N,1)*0.995);
for i = 1 : 99
    A(i,i+1) = 0.005;
end

H = Phi(:,:,1) * Psi;

Ce = zeros(N,N,nt);
Ce(:,:,1) = eye(N) * 10 ^ 10;
Cw = eye(N,N) * 1;
Cv = eye(nx) * 1;
ckv_k = Phi(:,:,1) * Cv * Phi(:,:,1)';

for k = 1 : nt-1
%     phi_k = Phi(:,:,k+1);
%     psi_k = phi_k * Psi;
%     ckv_k = phi_k * Cv * phi_k';
    psi_k = H;
    yk = y_sampling(:,k+1);
    zp = A * z(:,k);
    Cp = A * Ce(:,:,k) * A' + Cw;
    K = Cp * psi_k' * (psi_k * Cp * psi_k' + ckv_k)^(-1);
    z(:,k+1) = (eye(N) - K * psi_k) * zp + K * yk;
    Ce(:,:,k+1) = (1 - K * psi_k) * Cp * (1 - K * psi_k)';
end

y_re = zeros(nx,nt);
for i = 1 : nt
    y_re(:,i) = Psi * z(:,i);
end

for n = 1 : nt
    clf
    plot(x, y_re(:,n),'k-','LineWidth',1)
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





























