clc
clear
close all

load Messwerte.mat

M = zeros(8, 8 * 2001);

for i = 1 : 2001
    M(:, 8*(i-1)+1 : 8 * i) = reshape(m(:,i),8,8);
end

M = M / max(max(M));

imwrite(M,'Mess.png');