clc
clear
close all

printfigure = 0;

load('T_2D_noisefree.mat')

for i = 1 : 201
    f(:,i) = reshape(squeeze(F(i,:,:)),64*64,1);
end


f_fdm = zeros(64*64,201);
f_fdm(:,1) = f(:,1);

Dt = 0.1; % time_step

Lx = 10;
Ly = 20;

dx = 10 / 63;
x = 0 : dx : Lx;
nx = length(x);
dy = 20 / 63;
y = 0 : dy : Ly;
ny = length(y);

t = 0 : Dt : 20;
nt = length(t);

k = 0.1;  % Waermeleitfaehigkeit in cm^2/s

Knoten = @(m,n) (n-1) * nx + m;
Punkte = @(knoten) [mod(knoten,nx) (knoten-mod(knoten,nx))/nx+1];

dt = 0.1;
A = k * dt / dx^2;
B = k * dt / dy^2;
Q = 1 + 2*A + 2*B;

listQi = [];
listnAi = [];
listnBi = [];
listQAi = [];
listQBi = [];
listQABi = [];

listQj = [];
listnAj = [];
listnBj = [];
listQAj = [];
listQBj = [];
listQABj = [];


for m = 1 : nx
    for n = 1 : ny
        listQi = [listQi Knoten(m,n)];
        listnAi = [listnAi Knoten(m,n)];
        listnAi = [listnAi Knoten(m,n)];
        listnBi = [listnBi Knoten(m,n)];
        listnBi = [listnBi Knoten(m,n)];
        
        listQj = [listQj Knoten(m,n)];
        listnAj = [listnAj Knoten(min(m+1,nx),n)];
        listnAj = [listnAj Knoten(max(m-1,1),n)];
        listnBj = [listnBj Knoten(m,min(n+1,ny))];
        listnBj = [listnBj Knoten(m,max(n-1,1))];
    end
end


invD =   sparse(listQi,listQj,Q,nx*ny,nx*ny)...
    + sparse(listnAi,listnAj,-A,nx*ny,nx*ny)...
    + sparse(listnBi,listnBj,-B,nx*ny,nx*ny);

% D = inv(invD);
   

u = zeros(64*64, 201);

for i = 1 : 200
    u(:,i) = f(:,i+1) - invD \ f_fdm(:,i);
    f_fdm(:,i+1) = f(:,i+1);
end
    

for i = 1 : 201
    U(i,:,:) = reshape(u(:,i),64,64);
end

[Y,X] = meshgrid(y,x);

figure
for n = 1 : 5 : nt
    mesh(X,Y,squeeze(U(n,:,:)))
    zlim([-1 1])
    caxis([-1 3])
    pbaspect([1 Ly/Lx 0.5])
    setmesh('Temperature Distribution','$x$','$y$','$f$','T_2D_inhomo_modal_2',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(n)),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',30);
    set(TEXT,'Interpreter','latex')
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'T_2D_sr.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'T_2D_sr.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end

U = u;

save('Anregung.mat', 'U')









