function maketemperature(sigma_mu, sigma_sr, name)
printfigure = 0;

Length = 10;  % Stablaenge
Time = 20;   % Zetiraum

dx = Length / 1023;
dt = 0.01;
x = 0 : dx : Length;
t = 0 : dt : Time;
nx = length(x);
nt = length(t);

k = 0.1; % Waermeleitfaehigkeit in cm^2/s

N = 50;  % Grad
lambda = 0 : pi / Length : N * pi / Length;
f = zeros(nt, nx);  % Temperaturmatrix
f(1,:) = sin(x / Length * 2 * pi) + 1;
% f(1,:) = floor(x * 5) / 5;
phi = zeros(N + 1, nx);  % Eigenfunktionen
T = zeros(N + 1, nt);  % Gewichtung
u = zeros(nt, nx);  % Anregung
U = zeros(N + 1, nt);  % Anregungsgewichtung


phi(1,:) = sqrt(1 / Length);
for i = 2 : N + 1
    for n = 1 : nx
        phi(i, n) = sqrt(2 / Length) * cos(lambda(i) * x(n));
    end
end




u(:,round(3/dx)+1) = 0.1 * sin(t - pi / 4);
u(:,round(5/dx)+1) = -0.2 * sin(t);
u(:,round(7/dx)+1) = 0.01 * t;

for i = 1 : N + 1
    U(i,:) = u(:,round(3/dx)+1) * phi(i,round(3/dx)+1) + u(:,round(5/dx)+1) * phi(i,round(5/dx)+1) + u(:,round(7/dx)+1) * phi(i,round(7/dx)+1);
end

for i = 1 : N + 1
    T(i, 1) = 0;
    for n = 1 : nx
        T(i, 1) = T(i, 1) + f(1, n) * phi(i, n) * dx;
    end
end

A = zeros(N+1,N+1);
for i = 1 : N+1
    A(i,i) = (1 - dt * k * lambda(i)^2);
end


W = randn(N+1,nt-1) * sigma_sr;
V = randn(nx,nt) * sigma_mu;

for n = 2 : nt
    T(:, n) = A * T(:, n - 1) + dt * U(:, n-1) + W(:,n-1);
end

f = zeros(nt, nx);  % Temperaturmatrix
for n = 1 : nt
    for i = 1 : N + 1
        f(n,:) = f(n,:) + T(i,n) * phi(i,:);
    end
end

f_sr = f;
f_mu = f + V';


p_index = round(linspace(1,nx,64)); % 64 sensors from nx points
p = x(p_index);                     % positions of 51 sensors from 0 - Length
m = f_mu(:,p_index)';                  % their measurements

f_sr = f_sr';
f_mu = f_mu';

f_sr = f_sr(:,1:10:end);
f_mu = f_mu(:,1:10:end);

dt = 0.1;

m_rh = m(:,1:10:end);
save(name,'f_sr','f_mu','k','Length','dt','dx','p','m_rh','p_index','sigma_sr','sigma_mu')
end
