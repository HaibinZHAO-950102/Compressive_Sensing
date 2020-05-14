clc
clear
close all
printfigure = 0;

load('Messung')

f = y(:,end);
nt = length(t);
nx = length(x);

sampling_number = 40;
sampling_index = zeros(1,sampling_number);
sampling_index(1) = ceil(rand()*nx);
k = 2;
while k <= sampling_number
    temp = ceil(rand()*nx);
    if abs(sampling_index-temp) ~= 0
        sampling_index(k) = temp;
        k = k + 1;
    end
end
sampling_index = sort(sampling_index);

f_sampling = f(sampling_index);

Psi = zeros(nx,100);
for i = 1 : nx
    for j = 1 : 100
        Psi(i,j) = sin(2*pi*j*x(i));
    end
end

% Psi = idct(eye(nx));


Phi = zeros(sampling_number, nx);
for i = 1 : sampling_number
    Phi(i,sampling_index(i)) = 1;
end

A = Phi * Psi;

e = 0.001;
cvx_begin
    variable a(100,1)
    minimize(norm(a,1))
    subject to
        norm(A * a - f_sampling) <= e
cvx_end

f_re = Psi * a;
plot(x,f,'k-','LineWidth',1);
hold on
plot(x,f_re,'r-','LineWidth',1)
hold on
plot(x(sampling_index),f(sampling_index),'b.','markersize',10)
legend('Signal Distribution','Recovered Signal','Measurement')
% xlim([0 10])
% ylim([0 2.5])
setplt('Signal Recovery','$x$','$T$','Signal Recovery',printfigure)