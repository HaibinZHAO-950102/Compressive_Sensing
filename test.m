close all
printfigure = 1;

mesh(X,Y,squeeze(f(201,:,:)));
zlim([-1.5 1.5])
pbaspect([1 Length_y/Length_x 0.5])
set(gcf,'outerposition',get(0,'screensize'));
txt = ['$t = 200$'];
TEXT = text(8,0,0.6,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
txt = ['$N = ',num2str(N),'$'];
TEXT = text(8,0,1.5,txt,'FontSize',30);
set(TEXT,'Interpreter','latex')
setmesh('Temperature Distribution','$x$','$y$','$T$','T_2D_homo_modal_1_shot_5',1)