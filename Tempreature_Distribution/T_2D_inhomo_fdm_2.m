clc
clear
close all

printfigure = 1;

Length_x = 10;  
Length_y = 20;  
Time = 20;   % Zetiraum
dx = 0.25;
dy = 0.5;
dt = 0.4;
x = 0 : dx : Length_x;
y = 0 : dy : Length_y;
t = 0 : dt : Time;
nx = length(x);
ny = length(y);
nt = length(t);

k = 0.1;

f = zeros(nx, ny, nt);  % Temperaturmatrix
for m = 1 : nx
    for n = 1 : ny
        f(m,n,1) = sin(x(m) / Length_x * 2 * pi) * cos(y(n) / Length_y * 2 * pi) + 1;
    end
end

[Y,X] = meshgrid(y,x);
mesh(X,Y,squeeze(f(:,:,1)));
caxis([0 2.5])
pbaspect([1 Length_y/Length_x 0.5])
setmesh('Anfangsbedingung','$x$','$y$','$f$','T_2D_inhomo_fdm_inital_condition',1)

A = (dx ^ 2 * dy ^ 2) / (dx ^ 2 * dy ^ 2 + 2 * k * dt * (dx ^ 2 + dy ^ 2));
B = (k * dt * dx ^ 2) / (dx ^ 2 * dy ^ 2 + 2 * k * dt * (dx ^ 2 + dy ^ 2));
C = (k * dt * dy ^ 2) / (dx ^ 2 * dy ^ 2 + 2 * k * dt * (dx ^ 2 + dy ^ 2));
D = - (dt * dx ^ 2 * dy ^ 2) / (dx ^ 2 * dy ^ 2 + 2 * k * dt * (dx ^ 2 + dy ^ 2));

index = @(m,n,l) (l-1) * nx *ny + (n-1) * nx + m;    % Nummerierung

row = [];
colum = [];
value = [];
row_b = [];
colum_b = [];
value_b = [];

u = zeros(nx,ny,nt);
for l = 1 : nt
    u(3/dx+1,3/dy+1,l) = 1 * sin(t(l) - pi / 4) / dx / dy;
    u(5/dx+1,10/dy+1,l) = -2 * sin(t(l)) / dx / dy;
    u(7/dx+1,13/dy+1,l) = 0.02 * t(l) / dx / dy;
    u(8/dx+1,18/dy+1,l) = 0.05 * t(l) / dx / dy;
end


for m = 2 : nx - 1
    for n = 2 : ny - 1
        for l = 2 : nt
            m
            row = [row ; index(m,n,l); index(m,n,l); index(m,n,l); index(m,n,l); index(m,n,l); index(m,n,l)];
            colum = [colum; index(m,n,l); index(m,n,l-1); index(m,n-1,l); index(m,n+1,l); index(m-1,n,l); index(m+1,n,l)];
            value = [value; -1; A; B; B; C; C];
            row_b = [row_b; index(m,n,l)];
            colum_b = [colum_b; 1];
            value_b = [value_b; D * u(m,n,l)];
        end
    end
end

length_edge_x = [1, nx];
for m = 1 : 2
    for n = 1 : ny
        for l = 2 : nt
            m
            row = [row ; index(length_edge_x(m),n,l); index(length_edge_x(m),n,l)];
            colum = [colum; index(length_edge_x(m),n,l); index(length_edge_x(m)-sign(length_edge_x(m)-nx/2),n,l)];
            value = [value; -1; 1];
        end
    end
end

length_edge_y = [1, ny];
for m = 2 : nx - 1
    for n = 1 : 2
        for l = 2 : nt
            m
            row = [row ; index(m,length_edge_y(n),l); index(m,length_edge_y(n),l)];
            colum = [colum; index(m,length_edge_y(n),l); index(m,length_edge_y(n)-sign(length_edge_y(n)-ny/2),l)];
            value = [value; -1; 1];
        end
    end
end

for m = 1 : nx
    for n = 1 : ny
        m
        row = [row ; index(m,n,1)];
        colum = [colum; index(m,n,1)];
        value = [value; 1];
        row_b = [row_b; index(m,n,1)];
        colum_b = [colum_b; 1];
        value_b = [value_b; f(m,n,1)];
    end
end

H = sparse(row,colum,value,nx * ny * nt,nx * ny * nt);
b = sparse(row_b,colum_b,value_b,nx*ny*nt,1);

T = H\b;
for m = 1 : nx
    for n = 1 : ny
        for l = 2 : nt
            f(m,n,l) = T(index(m,n,l));
        end
    end
end

figure
for n = 1 : nt
    mesh(X,Y,squeeze(f(:,:,n)))
    zlim([-1 3])
    caxis([-1 3])
    pbaspect([1 Length_y/Length_x 0.5])
    setmesh('','$x$','$y$','$f$','T_2D_inhomo_fdm_1',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*dt),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    zticks([-1 1 3])
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    zticklabels({'-1','1','3'})
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'T_2D_inhomo_fdm_2.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'T_2D_inhomo_fdm_2.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end
    

figure
for n = 1 : 6
    a = round((nt - 1) / 5 * (n - 1)) + 1;
    mesh(X,Y,squeeze(f(:,:,a)))
    zlim([-1 3])
    caxis([-1 3])
    pbaspect([1 Length_y/Length_x 0.5])
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(a)),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    zticks([-1 1 3])
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    zticklabels({'-1','1','3'})
    drawnow
    name = ['T_2D_inhomo_fdm_2_shot_',num2str(n)];
    setmesh('','$x$','$y$','$f$',name,printfigure)
end

close all