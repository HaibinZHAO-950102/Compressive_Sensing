clc
clear
close all

printfigure = 1;

load('Messwerte_rh')

M = 24;  % Anzahl der Messungen
Dt = 0.1; % time_step

S = round(linspace(1,64,M));  % benutzte Sensoren
m_rh = m_rh(S,1:Dt/dt:end);


x = 0 : dx : 10;
nx = length(x);
t = 0 : Dt : 20;
nt = length(t);

k = 0.1;
E = sparse(2:nx,1:nx-1,1,nx,nx);
P = Dt*k/dx^2;
A = -2*speye(nx)+(E+E'); 
A(1,1) = -1; % Neumann Boundary
A(nx,nx) = -1; % Neumann Boundary
D = speye(nx) - P*A;


u = zeros(nt, nx);
u(:,round(3/dx+1)) = 0.1 * sin(t - pi / 4) / dx;
u(:,round(5/dx+1)) = -0.2 * sin(t) / dx;
u(:,round(7/dx+1)) = 0.01 * t / dx;




