clc
clear
close all

printfigure = 0;

I_origin = imread('Image.png');
I = imresize(I_origin,1/64);
imshow(I)

X = size(I,1);
Y = size(I,2);
N = X * Y;

I_vector = reshape(I, N, 1);

sampling_number = ceil(N * 0.5);

P_sampling = sort(ceil(rand(sampling_number, 1) * N));

I_sampling = I_vector(P_sampling);

I_sampling_2D = ones(N, 1) * 255;
I_sampling_2D(P_sampling) = I_vector(P_sampling);
I_sampling_2D = reshape(I_sampling_2D, X, Y);

figure
imshow(I_sampling_2D,[])
if printfigure == 1
    imwrite(I, 'Image_downsampling.png');
    imwrite(I_sampling_2D / 255, 'Image_Sampling.png');
end

psi_x = idct(eye(X));
psi_y = idct(eye(Y));

Psi = kron(psi_y,psi_x);
Phi = zeros(sampling_number, size(I,1) * size(I,2));
for i = 1 : sampling_number
    Phi(i,P_sampling(i)) = 1;
end
A = Phi * Psi;

cvx_begin
    variable a(N, 1)
    minimize(norm(a,1))
    subject to
        A * a - double(I_sampling) == zeros(size(A,1),1)
cvx_end

a_2D = reshape(a,X,Y);
I_re = idct2(a_2D);
I_re = I_re / max(max(I_re));
figure
imshow(I_re,[])
if printfigure == 1
    imwrite(I_re, 'Image_reconstruction.png');
end
     
                
                
                
                
                

