clc
clear

printfigure = 1;

Length = 10;  % Stablaenge
Time = 200;   % Zetiraum

step_length = 0.1;
step_time = 1;

x = 0 : step_length : Length;
t = 0 : step_time : Time;

N_length = length(x);
N_time = length(t);

k = 0.1;

f = zeros(N_time, N_length);
f(1,:) = sin(x / Length * 2 * pi);

A = step_length ^ 2 / (step_length ^ 2 + 2 * k * step_time);
B = (k * step_time) / (step_length ^ 2 + 2 * k * step_time);

M = speye(N_length * N_time);

for n = 2 : N_length - 1
    for m = 2 : N_time
        M(n+N_length*(m-1),n+N_length*(m-1)) = -1;
        M(n+N_length*(m-1),n+N_length*(m-2)) = A;
        M(n+N_length*(m-1),n+N_length*(m-1)-1) = B;
        M(n+N_length*(m-1),n+N_length*(m-1)+1) = B;
    end
end
length_edge = [1, N_length];
for n = 1 : 2
    for m = 2 : N_time
        M(length_edge(n)+N_length*(m-1),length_edge(n)+N_length*(m-1)) = -1;
        M(length_edge(n)+N_length*(m-1),length_edge(n)+N_length*(m-1)-sign(length_edge(n)-N_length/2)) = 1;
    end
end

b = zeros(N_length*N_time,1);
for n = 1 : N_length
    M(n,n) = 1;
    b(n) = f(1,n);
end

T = M\b;
for n = 1 : N_length
    for m = 2 : N_time
        f(m,n) = T(n+N_length*(m-1));
    end
end

figure
[X, Y] = meshgrid(x, t);
mesh(X,Y,f)
setmesh('Tempreature Distribution','$x$','$t$','$T$','TV_homo_fdm_Ttx_1',printfigure)
    