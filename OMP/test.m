clc
clear
close all

M = 12;
N = 64;
n = ceil(sqrt(N));
N = n^2;


trytime = 2500;
fehler = zeros(trytime,M);
rate = zeros(M,1);

for s = 1 : M
    s
for kkk = 1 : trytime


psi1 = idct(eye(n));
psi2 = idct(eye(n));

Psi = kron(psi1, psi2);

a = zeros(N,1);

sparsity = s;

a(randperm(N, sparsity)) = randn(sparsity,1) * 10;

y = Psi * a;

% plot(y,'k-','linewidth',5);

index = sort(randperm(N,M));
Phi = zeros(M,N);
for i = 1 : M
    Phi(i,index(i)) = 1;
end

H = Phi * Psi;

beo = Phi * y;

% hold on
% plot(index,beo,'b.','markersize',30)

x = OMP(H,beo,'STOPTOLERANCE',0.3);

% hold on
% plot(Psi * x,'r--','linewidth',5);
% xlim([0 N])
% setplt('Test OMP',' ','$Value$','test omp',0)

fehler(kkk,s) = norm(Psi*(x-a))/norm(Psi*a);
end
end

for sp = 1 : M
    a = find(fehler(:,sp) < 0.1);
    rate(sp) = length(a) / trytime;
end

rate

