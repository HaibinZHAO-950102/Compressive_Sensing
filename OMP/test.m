clc
clear
close all

M = 20;

t = 0 : 0.01 : 5;
x = sin(2*pi*2*t) + sin(2*pi*20*t);
plot(t,x)

f = 0 : 100;
for i = 1 : length(f)
    A(:,i) = sin(2*pi*f(i)*t);
end
A = [A,randn(501,400)];

sampling = randperm(length(t),M);

Phi = zeros(M,length(t));
for i = 1 : M
    Phi(i,sampling(i)) = 1;
end

a = OMP(2,Phi*x',Phi)