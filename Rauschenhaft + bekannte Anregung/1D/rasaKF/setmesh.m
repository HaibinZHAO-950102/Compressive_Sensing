function setmesh(figurename,xname,yname,zname,filename,printfigure)
    set(gca,'Fontsize',40)
    set(gca,'fontname','times new Roman')
    T = title(figurename,'fontsize',80);
    set(T,'Interpreter','latex')
    T = xlabel(xname,'fontsize',60);
    set(T,'Interpreter','latex')
    T = ylabel(yname,'fontsize',60);
    set(T,'Interpreter','latex')
    T = zlabel(zname,'fontsize',60);
    set(T,'Interpreter','latex')
    set(gcf,'outerposition',get(0,'screensize'));
    if printfigure == 1
        name = [filename,'.png'];
        exportgraphics(gcf, name,'Resolution',200)
    end
    drawnow
end