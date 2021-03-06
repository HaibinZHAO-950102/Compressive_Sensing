clc
clear
close all

printfigure = 1;

Length_x = 10;  
Length_y = 20;  
Time = 5;   % Zetiraum
step_length_x = 0.25;
step_length_y = 0.5;
step_time = 0.1;
x = 0 : step_length_x : Length_x;
y = 0 : step_length_y : Length_y;
t = 0 : step_time : Time;
N_length_x = length(x);
N_length_y = length(y);
N_time = length(t);

k = 0.1;
a = 0.1;

f = zeros(N_length_x, N_length_y, N_time);  % Temperaturmatrix
for m = 1 : N_length_x
    for n = 1 : N_length_y
        f(m,n,1) = sin(x(m) / Length_x * 2 * pi) * cos(y(n) / Length_y * 2 * pi) + 1;
    end
end

[Y,X] = meshgrid(y,x);
mesh(X,Y,squeeze(f(:,:,1)));
caxis([0 2.5])
pbaspect([1 Length_y/Length_x 0.5])
setmesh('Anfangsbedingung','$x$','$y$','$f$','T_2D_inhomo_fdm_inital_condition',printfigure)

A = (step_length_x ^ 2 * step_length_y ^ 2) / (step_length_x ^ 2 * step_length_y ^ 2 + 2 * k * step_time * (step_length_x ^ 2 + step_length_y ^ 2));
B = (k * step_time * step_length_x ^ 2) / (step_length_x ^ 2 * step_length_y ^ 2 + 2 * k * step_time * (step_length_x ^ 2 + step_length_y ^ 2));
C = (k * step_time * step_length_y ^ 2) / (step_length_x ^ 2 * step_length_y ^ 2 + 2 * k * step_time * (step_length_x ^ 2 + step_length_y ^ 2));
D = - (step_time * step_length_x ^ 2 * step_length_y ^ 2) / (step_length_x ^ 2 * step_length_y ^ 2 + 2 * k * step_time * (step_length_x ^ 2 + step_length_y ^ 2));
K = a * D - 1;

index = @(m,n,l) (l-1) * N_length_x *N_length_y + (n-1) * N_length_x + m;    % Nummerierung

row = [];
colum = [];
value = [];
row_b = [];
colum_b = [];
value_b = [];

for m = 2 : N_length_x - 1
    for n = 2 : N_length_y - 1
        for l = 2 : N_time
            row = [row ; index(m,n,l); index(m,n,l); index(m,n,l); index(m,n,l); index(m,n,l); index(m,n,l)];
            colum = [colum; index(m,n,l); index(m,n,l-1); index(m,n-1,l); index(m,n+1,l); index(m-1,n,l); index(m+1,n,l)];
            value = [value; K; A; B; B; C; C];
        end
    end
end

length_edge_x = [1, N_length_x];
for m = 1 : 2
    for n = 1 : N_length_y
        for l = 2 : N_time
            row = [row ; index(length_edge_x(m),n,l); index(length_edge_x(m),n,l)];
            colum = [colum; index(length_edge_x(m),n,l); index(length_edge_x(m)-sign(length_edge_x(m)-N_length_x/2),n,l)];
            value = [value; -1; 1];
        end
    end
end

length_edge_y = [1, N_length_y];
for m = 2 : N_length_x - 1
    for n = 1 : 2
        for l = 2 : N_time
            row = [row ; index(m,length_edge_y(n),l); index(m,length_edge_y(n),l)];
            colum = [colum; index(m,length_edge_y(n),l); index(m,length_edge_y(n)-sign(length_edge_y(n)-N_length_y/2),l)];
            value = [value; -1; 1];
        end
    end
end

for m = 1 : N_length_x
    for n = 1 : N_length_y
        row = [row ; index(m,n,1)];
        colum = [colum; index(m,n,1)];
        value = [value; 1];
        row_b = [row_b; index(m,n,1)];
        colum_b = [colum_b; 1];
        value_b = [value_b; f(m,n,1)];
    end
end

H = sparse(row,colum,value,N_length_x * N_length_y * N_time,N_length_x * N_length_y * N_time);
b = sparse(row_b,colum_b,value_b,N_length_x*N_length_y*N_time,1);

T = H\b;
for m = 1 : N_length_x
    for n = 1 : N_length_y
        for l = 2 : N_time
            f(m,n,l) = T(index(m,n,l));
        end
    end
end

figure
for n = 1 : N_time
    mesh(X,Y,squeeze(f(:,:,n)))
    zlim([0 2.5])
    caxis([0 2.5])
    pbaspect([1 Length_y/Length_x 0.5])
    setmesh('','$x$','$y$','$f$','T_2D_inhomo_fdm_1',0)
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str((n-1)*step_time),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    drawnow
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if printfigure == 1
        if n==1
             imwrite(imind,cm,'T_2D_inhomo_fdm_1.gif','gif', 'Loopcount',inf,'DelayTime',1e-6);
        else
             imwrite(imind,cm,'T_2D_inhomo_fdm_1.gif','gif','WriteMode','append','DelayTime',1e-6);
        end
    end
end
    

figure
for n = 1 : 6
    a = round((N_time - 1) / 5 * (n - 1)) + 1;
    mesh(X,Y,squeeze(f(:,:,a)))
    zlim([0 2.5])
    caxis([0 2.5])
    pbaspect([1 Length_y/Length_x 0.5])
    set(gcf,'outerposition',get(0,'screensize'));
    txt = ['$t = ',num2str(t(a)),'$'];
    TEXT = text(8,0,0.6,txt,'FontSize',60);
    set(TEXT,'Interpreter','latex')
    xticks(0:5:10)
    yticks(0:5:20)
    xticklabels({'0','5','10'})
    yticklabels({'0','5','10','15','20'})
    drawnow
    name = ['T_2D_inhomo_fdm_1_shot_',num2str(n)];
    setmesh('','$x$','$y$','$f$',name,printfigure)
end

close all