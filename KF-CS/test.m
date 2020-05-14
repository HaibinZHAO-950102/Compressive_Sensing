clc
clear
close all
printfigure = 0;

load('Messungen')

f = f_messung(end,:);
L = 10;
x = 0 : 0.01 : L;
N = length(x);

sampling_number = 100;
sampling_index = zeros(1,sampling_number);
sampling_index(1) = ceil(rand()*N);
k = 2;
while k <= sampling_number
    temp = ceil(rand()*N);
    if abs(sampling_index-temp) ~= 0
        sampling_index(k) = temp;
        k = k + 1;
    end
end
sampling_index = sort(sampling_index);

f_sampling = f(sampling_index);

% Psi = zeros(N, N);
% Psi(:,1) = sqrt(1/L);
% for i = 1 : N
%     for j = 2 : N
%         Psi(i,j) = sqrt(2/L)*cos(j*pi/L*x(i));
%     end
% end
% 
Psi = idct(eye(N));

Phi = zeros(sampling_number, N);
for i = 1 : sampling_number
    Phi(i,sampling_index(i)) = 1;
end

A = Phi * Psi;

e = 0.001;
cvx_begin
    variable a(N,1)
    minimize(norm(a,1))
    subject to
        norm(A * a - f_sampling') <= e
cvx_end

f_re = Psi * a;
plot(x,f,'k-','LineWidth',1);
hold on
plot(x,f_re,'r-','LineWidth',1)
hold on
plot(x(sampling_index),f(sampling_index),'b.','markersize',10)
legend('Tempreature Distribution','Recovered Tempreature','Measurement')
txt = ['$N = ',num2str(sampling_number),'$'];
T = text(0.8,0.8,txt,'FontSize',30);
set(T,'Interpreter','latex')
xlim([0 10])
ylim([0 2.5])
setplt('Tempreature Recovery','$x$','$T$','Tempreature Recovery',printfigure)
