% CS Basis
clc
clear

load Messwerte


DCT = idct(eye(64))^-1;



WT = DWT(64,'haar');



NT = WT * randn(64);
NT = orth(NT);



lambda = [0 : 63] * pi / 10;
Psi = zeros(64,64);
Psi(:, 1) = sqrt(1 / 10);
for i = 2 : 64
    for n = 1 : 64
        Psi(n, i) = sqrt(2 / 10) * cos(lambda(i) * p(n));
    end
end
Psi = orth(Psi);
CT = Psi^-1;

FFT = ifft(eye(64))^-1;

MT = WT*FFT;
THETA = MT^-1;

save('basis.mat','MT','THETA')