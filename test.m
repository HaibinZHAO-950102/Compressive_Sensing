clc
clear
load('Messwerte.mat')

% 'db1' or 'haar', 'db2', ... ,'db10', ... , 'db45'
% 'coif1', ... , 'coif5'
% 'sym2', ... , 'sym8', ... ,'sym45'
% 'bior1.1', 'bior1.3', 'bior1.5'
% 'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8'
% 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7'
% 'bior3.9', 'bior4.4', 'bior5.5', 'bior6.8'
% 'rbio1.1', 'rbio1.3', 'rbio1.5'
% 'rbio2.2', 'rbio2.4', 'rbio2.6', 'rbio2.8'
% 'rbio3.1', 'rbio3.3', 'rbio3.5', 'rbio3.7'
% 'rbio3.9', 'rbio4.4', 'rbio5.5', 'rbio6.8'

y = f(end,1:end-1);
x = 0 : 0.01 : 10-0.001;
N = length(y);

Psi = DWT(N, 'haar')^-1;

M = 100;
sampling_index = sort(ceil(rand(1,M)*N));

Phi = zeros(M,N);
for i = 1 : M
    Phi(i,sampling_index(i)) = 1;
end

A = Phi * Psi;
y_sampling = y(sampling_index);

e = 0.01;
cvx_begin
    variable a(N,1)
    minimize(norm(a,1))
    subject to
        norm(A * a - y_sampling') <= e
cvx_end

y_re = Psi * a;
plot(x,y)
hold on
plot(x,y_re)
hold on
plot(x(sampling_index),y_sampling,'*')