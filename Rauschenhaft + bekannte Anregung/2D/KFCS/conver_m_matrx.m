clc
clear
close all

load Messwerte_rh

M = zeros(8, 8 * 201);

for i = 1 : 201
    M(:,8*(i-1)+1:8*i) = reshape(m_rh(:,i),8,8);
end

M = M / max(max(M));

imwrite(M,'Mess2D.png')

