clc
clear
close all

load OSB

scale = 100;

a = max(max(THETA));
b = min(min(THETA));

THETA = (THETA - b) / (a-b);

show = zeros(size(THETA) * scale);

for i = 1 : size(THETA,1)
    for j = 1 : size(THETA,2)
        show((i-1)*scale+1:i*scale,(j-1)*scale+1:j*scale) = THETA(i,j);
    end
end

figure
imshow(show)
imwrite(show,'OSB.png')

koe = zeros(size(THETA,2),1);
index = randperm(64,3);
koe(index) = abs(randn(3,1));

a = max(max(koe));
b = min(min(koe));

koe = (koe - b) / (a-b);
scale = 200;

show_koe = zeros(size(koe) * scale);

for i = 1 : size(koe,1)
    for j = 1 : size(koe,2)
        show_koe((i-1)*scale+1:i*scale,(j-1)*scale+1:j*scale) = koe(i,j);
    end
end

figure
imshow(show_koe)
imwrite(show_koe,'show_coefficient.png')


index = randperm(64,30);
phi = zeros(30,64);
for i = 1 : 30
    phi(i,index(i)) = 1;
end

THETA_phi = phi * THETA;
y = phi * THETA * koe;

    
a = max(max(THETA_phi));
b = min(min(THETA_phi));

THETA_phi = (THETA_phi - b) / (a-b);
scale = 100;

show_theta = zeros(size(THETA_phi) * scale);

for i = 1 : size(THETA_phi,1)
    for j = 1 : size(THETA_phi,2)
        show_theta((i-1)*scale+1:i*scale,(j-1)*scale+1:j*scale) = THETA_phi(i,j);
    end
end

figure
imshow(show_theta)
imwrite(show_theta,'show_theta.png')

    
    
a = max(max(y));
b = min(min(y));

y = (y - b) / (a-b);
scale = 200;

show_y = zeros(size(y) * scale);

for i = 1 : size(y,1)
    for j = 1 : size(y,2)
        show_y((i-1)*scale+1:i*scale,(j-1)*scale+1:j*scale) = y(i,j);
    end
end

figure
imshow(show_y)
imwrite(show_y,'show_y.png')

    
    
    
    
    
    
    
    
    
    
    
    
    
    

